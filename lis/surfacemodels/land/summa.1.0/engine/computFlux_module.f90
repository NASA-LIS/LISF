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

module computFlux_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! layer types
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! access the global print flag
USE globalData,only:globalPrintFlag

! control parameters
USE globalData,only:verySmall       ! a very small number
USE globalData,only:veryBig         ! a very big number
USE globalData,only:dx              ! finite difference increment

! constants
USE multiconst,only:&
                    gravity,      & ! acceleration of gravity              (m s-2)
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:       &
 localColumn,                     & ! separate groundwater representation in each local soil column
 singleBasin                        ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
 qbaseTopmodel,                   & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                       & ! a big bucket (lumped aquifer model)
 noExplicit                         ! no explicit groundwater parameterization

! look-up values for the form of Richards' equation
USE mDecisions_module,only:       &
 moisture,                        & ! moisture-based form of Richards' equation
 mixdform                           ! mixed form of Richards' equation

! look-up values for the choice of boundary conditions for hydrology
USE mDecisions_module,only:       &
 prescribedHead,                  & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 funcBottomHead,                  & ! function of matric head in the lower-most layer
 freeDrainage,                    & ! free drainage
 liquidFlux,                      & ! liquid water flux
 zeroFlux                           ! zero flux

implicit none
private
public::computFlux
public::soilCmpres
contains

 ! *********************************************************************************************************
 ! public subroutine computFlux: compute model fluxes
 ! *********************************************************************************************************
 subroutine computFlux(&
                       ! input-output: model control
                       nSnow,                   & ! intent(in):    number of snow layers
                       nSoil,                   & ! intent(in):    number of soil layers
                       nLayers,                 & ! intent(in):    total number of layers
                       firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                       firstFluxCall,           & ! intent(inout): flag to denote the first flux call
                       firstSplitOper,          & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                       computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                       drainageMeltPond,        & ! intent(in):    drainage from the surface melt pond (kg m-2 s-1)
                       ! input: state variables
                       scalarCanairTempTrial,   & ! intent(in):    trial value for the temperature of the canopy air space (K)
                       scalarCanopyTempTrial,   & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                       mLayerTempTrial,         & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                       mLayerMatricHeadLiqTrial,& ! intent(in):    trial value for the liquid water matric potential in each soil layer (m)
                       ! input: diagnostic variables defining the liquid water and ice content
                       scalarCanopyLiqTrial,    & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                       scalarCanopyIceTrial,    & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                       mLayerVolFracLiqTrial,   & ! intent(in):    trial value for the volumetric liquid water content in each snow and soil layer (-)
                       mLayerVolFracIceTrial,   & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                       ! input: data structures
                       model_decisions,         & ! intent(in):    model decisions
                       type_data,               & ! intent(in):    type of vegetation and soil
                       attr_data,               & ! intent(in):    spatial attributes
                       mpar_data,               & ! intent(in):    model parameters
                       forc_data,               & ! intent(in):    model forcing data
                       bvar_data,               & ! intent(in):    average model variables for the entire basin
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       indx_data,               & ! intent(in):    index data
                       ! input-output: data structures
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,               & ! intent(inout): model fluxes for a local HRU
                       deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ! input-output: flux vector and baseflow derivatives
                       ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                       fluxVec,                 & ! intent(out):   flux vector (mixed units)
                       ! output: error control
                       err,message)               ! intent(out):   error code and error message
 ! provide access to soil utilities
 !USE snow_utils_module,only:dFracLiq_dTk              ! differentiate the freezing curve w.r.t. temperature (snow)
 !USE soil_utils_module,only:dTheta_dPsi               ! derivative in the soil water characteristic (soil)
 !USE soil_utils_module,only:dPsi_dTheta               ! derivative in the soil water characteristic (soil)
 !USE soil_utils_module,only:dTheta_dTk                ! differentiate the freezing curve w.r.t. temperature (soil)
 !USE soil_utils_module,only:matricHead                ! compute the matric head based on volumetric water content
 ! provide access to flux subroutines
 USE vegNrgFlux_module,only:vegNrgFlux            ! compute energy fluxes over vegetation
 USE ssdNrgFlux_module,only:ssdNrgFlux            ! compute energy fluxes throughout the snow and soil subdomains
 USE vegLiqFlux_module,only:vegLiqFlux            ! compute liquid water fluxes through vegetation
 USE snowLiqFlx_module,only:snowLiqflx            ! compute liquid water fluxes through snow
 USE soilLiqFlx_module,only:soilLiqflx            ! compute liquid water fluxes through soil
 USE groundwatr_module,only:groundwatr            ! compute the baseflow flux
 ! provide access to indices that define elements of the data structures
 USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
 USE var_lookup,only:iLookPARAM                   ! named variables for structure elements
 USE var_lookup,only:iLookFORCE                   ! named variables for structure elements
 USE var_lookup,only:iLookPROG                    ! named variables for structure elements
 USE var_lookup,only:iLookINDEX                   ! named variables for structure elements
 USE var_lookup,only:iLookDIAG                    ! named variables for structure elements
 USE var_lookup,only:iLookFLUX                    ! named variables for structure elements
 USE var_lookup,only:iLookDERIV                   ! named variables for structure elements
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input-output: control
 integer(i4b),intent(in)         :: nSnow                       ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                       ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                     ! total number of layers
 logical(lgt),intent(in)         :: firstSubStep                ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall               ! flag to indicate if we are processing the first flux call
 logical(lgt),intent(in)         :: firstSplitOper              ! flag to indicate if we are processing the first flux call in a splitting operation
 logical(lgt),intent(in)         :: computeVegFlux              ! flag to indicate if computing fluxes over vegetation
 real(dp),intent(in)             :: drainageMeltPond            ! drainage from the surface melt pond (kg m-2 s-1)
 ! input: state variables
 real(dp),intent(in)             :: scalarCanairTempTrial       ! trial value for temperature of the canopy air space (K)
 real(dp),intent(in)             :: scalarCanopyTempTrial       ! trial value for temperature of the vegetation canopy (K)
 real(dp),intent(in)             :: mLayerTempTrial(:)          ! trial value for temperature of each snow/soil layer (K)
 real(dp),intent(in)             :: mLayerMatricHeadLiqTrial(:) ! trial value for the liquid water matric potential (m)
 ! input: diagnostic variables
 real(dp),intent(in)             :: scalarCanopyLiqTrial        ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)             :: scalarCanopyIceTrial        ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)             :: mLayerVolFracLiqTrial(:)    ! trial value for volumetric fraction of liquid water (-)
 real(dp),intent(in)             :: mLayerVolFracIceTrial(:)    ! trial value for volumetric fraction of ice (-)
 ! input: data structures
 type(model_options),intent(in)  :: model_decisions(:)          ! model decisions
 type(var_i),        intent(in)  :: type_data                   ! type of vegetation and soil
 type(var_d),        intent(in)  :: attr_data                   ! spatial attributes
 type(var_dlength),  intent(in)  :: mpar_data                   ! model parameters
 type(var_d),        intent(in)  :: forc_data                   ! model forcing data
 type(var_dlength),  intent(in)  :: bvar_data                   ! model variables for the local basin
 type(var_dlength),  intent(in)  :: prog_data                   ! prognostic variables for a local HRU
 type(var_ilength),  intent(in)  :: indx_data                   ! indices defining model states and layers
 ! input-output: data structures
 type(var_dlength),intent(inout) :: diag_data                   ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                   ! model fluxes for a local HRU
 type(var_dlength),intent(inout) :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
 ! input-output: flux vector and baseflow derivatives
 integer(i4b),intent(inout)      :: ixSaturation                ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(dp),intent(out)            :: dBaseflow_dMatric(:,:)      ! derivative in baseflow w.r.t. matric head (s-1)
 real(dp),intent(out)            :: fluxVec(:)                  ! model flux vector (mixed units)
 ! output: error control
 integer(i4b),intent(out)        :: err                         ! error code
 character(*),intent(out)        :: message                     ! error message
 ! ---------------------------------------------------------------------------------------
 ! * local variables
 ! ---------------------------------------------------------------------------------------
 integer(i4b)                    :: local_ixGroundwater         ! local index for groundwater representation
 integer(i4b)                    :: iLayer                      ! index of model layers
 real(dp),dimension(nSoil)       :: dHydCond_dMatric            ! derivative in hydraulic conductivity w.r.t matric head (s-1)
 character(LEN=256)              :: cmessage                    ! error message of downwind routine
 ! --------------------------------------------------------------
 ! initialize error control
 err=0; message='computFlux/'

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! get the necessary variables for the flux computations
 associate(&

 ! model decisions
 ixGroundwater                => model_decisions(iLookDECISIONS%groundwatr)%iDecision            ,& ! intent(in): [i4b]    groundwater parameterization
 ixSpatialGroundwater         => model_decisions(iLookDECISIONS%spatial_gw)%iDecision            ,& ! intent(in): [i4b]    spatial representation of groundwater (local-column or single-basin)

 ! domain boundary conditions
 upperBoundTemp               => forc_data%var(iLookFORCE%airtemp)                               ,& ! intent(in): [dp]     temperature of the upper boundary of the snow and soil domains (K)
 scalarRainfall               => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)                  ,& ! intent(in): [dp]     rainfall rate (kg m-2 s-1)

 ! canopy and layer depth
 canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)               ,& ! intent(in): [dp   ]  canopy depth (m)
 mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat                        ,& ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)

 ! indices of model state variables for the vegetation subdomain
 ixCasNrg                     => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)                       ,& ! intent(in): [i4b]    index of canopy air space energy state variable
 ixVegNrg                     => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)                       ,& ! intent(in): [i4b]    index of canopy energy state variable
 ixVegHyd                     => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)                       ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
 ixTopNrg                     => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)                       ,& ! intent(in): [i4b]    index of upper-most energy state in the snow+soil subdomain
 ixTopHyd                     => indx_data%var(iLookINDEX%ixTopHyd)%dat(1)                       ,& ! intent(in): [i4b]    index of upper-most hydrology state in the snow+soil subdomain

 ! indices of model state variables for the snow+soil domain
 ixSnowSoilNrg                => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat                     ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
 ixSnowSoilHyd                => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat                     ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
 layerType                    => indx_data%var(iLookINDEX%layerType)%dat                        , & ! intent(in): [i4b(:)] type of layer (iname_soil or iname_snow)

 ! number of state variables of a specific type
 nSnowSoilNrg                 => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
 nSnowOnlyNrg                 => indx_data%var(iLookINDEX%nSnowOnlyNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
 nSoilOnlyNrg                 => indx_data%var(iLookINDEX%nSoilOnlyNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
 nSnowSoilHyd                 => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
 nSnowOnlyHyd                 => indx_data%var(iLookINDEX%nSnowOnlyHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
 nSoilOnlyHyd                 => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain

 ! snow parameters
 snowfrz_scale                => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)                  ,& ! intent(in): [dp] scaling parameter for the snow freezing curve (K-1)

 ! derivatives
 dPsiLiq_dPsi0                => deriv_data%var(iLookDERIV%dPsiLiq_dPsi0   )%dat                 ,&  ! intent(in):  [dp(:)] derivative in liquid water matric pot w.r.t. the total water matric pot (-)
 dPsiLiq_dTemp                => deriv_data%var(iLookDERIV%dPsiLiq_dTemp   )%dat                 ,&  ! intent(in):  [dp(:)] derivative in the liquid water matric potential w.r.t. temperature
 mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat                 ,&  ! intent(in):  [dp(:)] derivative of volumetric liquid water content w.r.t. temperature
 dTheta_dTkCanopy             => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)              ,&  ! intent(in):  [dp]    derivative of volumetric liquid water content w.r.t. temperature

 ! number of flux calls
 numFluxCalls                 => diag_data%var(iLookDIAG%numFluxCalls)%dat(1)                    ,&  ! intent(out): [dp] number of flux calls (-) 

 ! net fluxes over the vegetation domain
 scalarCanairNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1)          ,&  ! intent(out): [dp] net energy flux for the canopy air space        (W m-2)
 scalarCanopyNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1)          ,&  ! intent(out): [dp] net energy flux for the vegetation canopy       (W m-2)
 scalarGroundNetNrgFlux       => flux_data%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1)          ,&  ! intent(out): [dp] net energy flux for the ground surface          (W m-2)
 scalarCanopyNetLiqFlux       => flux_data%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1)          ,&  ! intent(out): [dp] net liquid water flux for the vegetation canopy (kg m-2 s-1)

 ! net fluxes over the snow+soil domain
 mLayerNrgFlux                => flux_data%var(iLookFLUX%mLayerNrgFlux)%dat                      ,&  ! intent(out): [dp] net energy flux for each layer within the snow+soil domain (J m-3 s-1)
 mLayerLiqFluxSnow            => flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat                  ,&  ! intent(out): [dp] net liquid water flux for each snow layer (s-1)
 mLayerLiqFluxSoil            => flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat                  ,&  ! intent(out): [dp] net liquid water flux for each soil layer (s-1)

 ! evaporative fluxes
 scalarCanopyTranspiration    => flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1)       ,&  ! intent(out): [dp]    canopy transpiration (kg m-2 s-1)
 scalarCanopyEvaporation      => flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)         ,&  ! intent(out): [dp]    canopy evaporation/condensation (kg m-2 s-1)
 scalarGroundEvaporation      => flux_data%var(iLookFLUX%scalarGroundEvaporation)%dat(1)         ,&  ! intent(out): [dp]    ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 mLayerTranspire              => flux_data%var(iLookFLUX%mLayerTranspire)%dat                    ,&  ! intent(out): [dp(:)] transpiration loss from each soil layer (m s-1)

 ! fluxes for the snow+soil domain
 iLayerNrgFlux                => flux_data%var(iLookFLUX%iLayerNrgFlux)%dat                      ,&  ! intent(out): [dp(0:)] vertical energy flux at the interface of snow and soil layers
 iLayerLiqFluxSnow            => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat                  ,&  ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
 iLayerLiqFluxSoil            => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat                  ,&  ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
 mLayerHydCond                => flux_data%var(iLookFLUX%mLayerHydCond)%dat                      ,&  ! intent(out): [dp(:)]  hydraulic conductivity in each soil layer (m s-1)
 mLayerBaseflow               => flux_data%var(iLookFLUX%mLayerBaseflow)%dat                     ,&  ! intent(out): [dp(:)]  baseflow from each soil layer (m s-1)
 scalarSnowDrainage           => flux_data%var(iLookFLUX%scalarSnowDrainage)%dat(1)              ,&  ! intent(out): [dp]     drainage from the snow profile (m s-1)
 scalarSoilDrainage           => flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1)              ,&  ! intent(out): [dp]     drainage from the soil profile (m s-1)
 scalarSoilBaseflow           => flux_data%var(iLookFLUX%scalarSoilBaseflow)%dat(1)              ,&  ! intent(out): [dp]     total baseflow from the soil profile (m s-1)

 ! infiltration
 scalarInfilArea              => diag_data%var(iLookDIAG%scalarInfilArea   )%dat(1)              ,&  ! intent(out): [dp] fraction of unfrozen area where water can infiltrate (-)
 scalarFrozenArea             => diag_data%var(iLookDIAG%scalarFrozenArea  )%dat(1)              ,&  ! intent(out): [dp] fraction of area that is considered impermeable due to soil ice (-)
 scalarSoilControl            => diag_data%var(iLookDIAG%scalarSoilControl )%dat(1)              ,&  ! intent(out): [dp] soil control on infiltration, zero or one
 scalarMaxInfilRate           => flux_data%var(iLookFLUX%scalarMaxInfilRate)%dat(1)              ,&  ! intent(out): [dp] maximum infiltration rate (m s-1)
 scalarInfiltration           => flux_data%var(iLookFLUX%scalarInfiltration)%dat(1)              ,&  ! intent(out): [dp] infiltration of water into the soil profile (m s-1)

 ! boundary fluxes in the soil domain
 scalarThroughfallRain        => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)           ,&  ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 scalarCanopyLiqDrainage      => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)         ,&  ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 scalarRainPlusMelt           => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1)              ,&  ! intent(out): [dp] rain plus melt (m s-1)
 scalarSurfaceRunoff          => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)             ,&  ! intent(out): [dp] surface runoff (m s-1)
 scalarExfiltration           => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1)              ,&  ! intent(out): [dp] exfiltration from the soil profile (m s-1)
 mLayerColumnOutflow          => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat                ,&  ! intent(out): [dp(:)] column outflow from each soil layer (m3 s-1)

 ! fluxes for the aquifer
 scalarAquiferTranspire       => flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1)          ,&  ! intent(out): [dp] transpiration loss from the aquifer (m s-1
 scalarAquiferRecharge        => flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1)           ,&  ! intent(out): [dp] recharge to the aquifer (m s-1)
 scalarAquiferBaseflow        => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)           ,&  ! intent(out): [dp] total baseflow from the aquifer (m s-1)

 ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
 dCanairNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanairTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy air space flux w.r.t. canopy air temperature
 dCanairNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanopyTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy air space flux w.r.t. canopy temperature
 dCanairNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dGroundTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy air space flux w.r.t. ground temperature
 dCanopyNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanairTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy flux w.r.t. canopy air temperature
 dCanopyNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanopyTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy flux w.r.t. canopy temperature
 dCanopyNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dGroundTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy flux w.r.t. ground temperature
 dCanopyNetFlux_dCanLiq       => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanLiq      )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy fluxes w.r.t. canopy liquid water content
 dGroundNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanairTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net ground flux w.r.t. canopy air temperature
 dGroundNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanopyTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net ground flux w.r.t. canopy temperature
 dGroundNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dGroundTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net ground flux w.r.t. ground temperature
 dGroundNetFlux_dCanLiq       => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanLiq      )%dat(1)  ,&  ! intent(out): [dp] derivative in net ground fluxes w.r.t. canopy liquid water content

 ! derivatives in evaporative fluxes w.r.t. relevant state variables
 dCanopyEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanair )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy evaporation w.r.t. canopy air temperature
 dCanopyEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanopy )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy evaporation w.r.t. canopy temperature
 dCanopyEvaporation_dTGround  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTGround )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy evaporation w.r.t. ground temperature
 dCanopyEvaporation_dCanLiq   => deriv_data%var(iLookDERIV%dCanopyEvaporation_dCanLiq  )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy evaporation w.r.t. canopy liquid water content
 dGroundEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanair )%dat(1)  ,&  ! intent(out): [dp] derivative in ground evaporation w.r.t. canopy air temperature
 dGroundEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanopy )%dat(1)  ,&  ! intent(out): [dp] derivative in ground evaporation w.r.t. canopy temperature
 dGroundEvaporation_dTGround  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTGround )%dat(1)  ,&  ! intent(out): [dp] derivative in ground evaporation w.r.t. ground temperature
 dGroundEvaporation_dCanLiq   => deriv_data%var(iLookDERIV%dGroundEvaporation_dCanLiq  )%dat(1)  ,&  ! intent(out): [dp] derivative in ground evaporation w.r.t. canopy liquid water content

 ! derivatives in canopy water w.r.t canopy temperature
 dCanLiq_dTcanopy             => deriv_data%var(iLookDERIV%dCanLiq_dTcanopy            )%dat(1)  ,&  ! intent(out): [dp] derivative of canopy liquid storage w.r.t. temperature

 ! derivatives in canopy liquid fluxes w.r.t. canopy water
 scalarCanopyLiqDeriv         => deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv        )%dat(1)  ,&  ! intent(out): [dp] derivative in (throughfall + drainage) w.r.t. canopy liquid water
 scalarThroughfallRainDeriv   => deriv_data%var(iLookDERIV%scalarThroughfallRainDeriv  )%dat(1)  ,&  ! intent(out): [dp] derivative in throughfall w.r.t. canopy liquid water
 scalarCanopyLiqDrainageDeriv => deriv_data%var(iLookDERIV%scalarCanopyLiqDrainageDeriv)%dat(1)  ,&  ! intent(out): [dp] derivative in canopy drainage w.r.t. canopy liquid water

 ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
 dNrgFlux_dTempAbove          => deriv_data%var(iLookDERIV%dNrgFlux_dTempAbove         )%dat     ,&  ! intent(out): [dp(:)] derivatives in the flux w.r.t. temperature in the layer above
 dNrgFlux_dTempBelow          => deriv_data%var(iLookDERIV%dNrgFlux_dTempBelow         )%dat     ,&  ! intent(out): [dp(:)] derivatives in the flux w.r.t. temperature in the layer below

 ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
 iLayerLiqFluxSnowDeriv       => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv      )%dat     ,&  ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces

 ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
 dVolTot_dPsi0                => deriv_data%var(iLookDERIV%dVolTot_dPsi0               )%dat     ,&  ! intent(out): [dp(:)] derivative in total water content w.r.t. total water matric potential
 dq_dHydStateAbove            => deriv_data%var(iLookDERIV%dq_dHydStateAbove           )%dat     ,&  ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
 dq_dHydStateBelow            => deriv_data%var(iLookDERIV%dq_dHydStateBelow           )%dat     ,&  ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
 mLayerdTheta_dPsi            => deriv_data%var(iLookDERIV%mLayerdTheta_dPsi           )%dat     ,&  ! intent(out): [dp(:)] derivative in the soil water characteristic w.r.t. psi
 mLayerdPsi_dTheta            => deriv_data%var(iLookDERIV%mLayerdPsi_dTheta           )%dat     ,&  ! intent(out): [dp(:)] derivative in the soil water characteristic w.r.t. theta
 dCompress_dPsi               => deriv_data%var(iLookDERIV%dCompress_dPsi              )%dat     ,&  ! intent(out): [dp(:)] derivative in compressibility w.r.t matric head

 ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
 dq_dNrgStateAbove            => deriv_data%var(iLookDERIV%dq_dNrgStateAbove           )%dat     ,&  ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
 dq_dNrgStateBelow            => deriv_data%var(iLookDERIV%dq_dNrgStateBelow           )%dat      &  ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
 )  ! association to data in structures

 ! *****
 ! * PRELIMINARIES... 
 ! ******************

 ! increment the number of flux calls
 numFluxCalls = numFluxCalls+1

 ! modify the groundwater representation for this single-column implementation
 select case(ixSpatialGroundwater)
  case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation)

 ! initialize liquid water fluxes throughout the snow and soil domains
 ! NOTE: used in the energy routines, which is called before the hydrology routines
 if(firstFluxCall)then
  if(nSnow > 0) iLayerLiqFluxSnow(0:nSnow) = 0._dp
                iLayerLiqFluxSoil(0:nSoil) = 0._dp
 end if

 ! *****
 ! * CALCULATE ENERGY FLUXES OVER VEGETATION...
 ! *********************************************

 ! check if there is a need to calculate the energy fluxes over vegetation
 if(ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing .or. ixTopNrg/=integerMissing .or. firstFluxCall)then

  ! derivative in canopy liquid storage w.r.t. canopy temperature
  dCanLiq_dTcanopy = dTheta_dTkCanopy*iden_water*canopyDepth  ! kg m-2 K-1

  ! calculate the energy fluxes over vegetation
  call vegNrgFlux(&
                  ! input: model control
                  firstSubStep,                           & ! intent(in): flag to indicate if we are processing the first sub-step
                  firstFluxCall,                          & ! intent(in): flag to indicate if we are processing the first flux call
                  computeVegFlux,                         & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                  ! input: model state variables
                  upperBoundTemp,                         & ! intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
                  scalarCanairTempTrial,                  & ! intent(in): trial value of the canopy air space temperature (K)
                  scalarCanopyTempTrial,                  & ! intent(in): trial value of canopy temperature (K)
                  mLayerTempTrial(1),                     & ! intent(in): trial value of ground temperature (K)
                  scalarCanopyIceTrial,                   & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                  scalarCanopyLiqTrial,                   & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                  ! input: model derivatives
                  dCanLiq_dTcanopy,                       & ! intent(in): derivative in canopy liquid storage w.r.t. canopy temperature (kg m-2 K-1)
                  ! input/output: data structures
                  type_data,                              & ! intent(in):    type of vegetation and soil
                  attr_data,                              & ! intent(in):    spatial attributes
                  forc_data,                              & ! intent(in):    model forcing data
                  mpar_data,                              & ! intent(in):    model parameters
                  indx_data,                              & ! intent(in):    index data
                  prog_data,                              & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                              & ! intent(inout): model diagnostic variables for a local HRU
                  flux_data,                              & ! intent(inout): model fluxes for a local HRU
                  bvar_data,                              & ! intent(in):    model variables for the local basin
                  model_decisions,                        & ! intent(in):    model decisions
                  ! output: liquid water fluxes associated with evaporation/transpiration
                  scalarCanopyTranspiration,              & ! intent(out): canopy transpiration (kg m-2 s-1)
                  scalarCanopyEvaporation,                & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                  scalarGroundEvaporation,                & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                  ! output: fluxes
                  scalarCanairNetNrgFlux,                 & ! intent(out): net energy flux for the canopy air space (W m-2)
                  scalarCanopyNetNrgFlux,                 & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                  scalarGroundNetNrgFlux,                 & ! intent(out): net energy flux for the ground surface (W m-2)
                  ! output: flux derivatives
                  dCanairNetFlux_dCanairTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                  dCanairNetFlux_dCanopyTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                  dCanairNetFlux_dGroundTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                  dCanopyNetFlux_dCanairTemp,             & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                  dCanopyNetFlux_dCanopyTemp,             & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                  dCanopyNetFlux_dGroundTemp,             & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                  dGroundNetFlux_dCanairTemp,             & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                  dGroundNetFlux_dCanopyTemp,             & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                  dGroundNetFlux_dGroundTemp,             & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                  ! output: liquid water flux derivarives (canopy evap)
                  dCanopyEvaporation_dCanLiq,             & ! intent(out): derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
                  dCanopyEvaporation_dTCanair,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                  dCanopyEvaporation_dTCanopy,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
                  dCanopyEvaporation_dTGround,            & ! intent(out): derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
                  ! output: liquid water flux derivarives (ground evap)
                  dGroundEvaporation_dCanLiq,             & ! intent(out): derivative in ground evaporation w.r.t. canopy liquid water content (s-1)
                  dGroundEvaporation_dTCanair,            & ! intent(out): derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                  dGroundEvaporation_dTCanopy,            & ! intent(out): derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
                  dGroundEvaporation_dTGround,            & ! intent(out): derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
                  ! output: cross derivative terms
                  dCanopyNetFlux_dCanLiq,                 & ! intent(out): derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                  dGroundNetFlux_dCanLiq,                 & ! intent(out): derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! check fluxes
  if(globalPrintFlag)then
   print*, '**'
   write(*,'(a,1x,10(f30.20))') 'canopyDepth           = ',  canopyDepth
   write(*,'(a,1x,10(f30.20))') 'mLayerDepth(1:2)      = ',  mLayerDepth(1:2)
   write(*,'(a,1x,10(f30.20))') 'scalarCanairTempTrial = ',  scalarCanairTempTrial   ! trial value of the canopy air space temperature (K)
   write(*,'(a,1x,10(f30.20))') 'scalarCanopyTempTrial = ',  scalarCanopyTempTrial   ! trial value of canopy temperature (K)
   write(*,'(a,1x,10(f30.20))') 'mLayerTempTrial(1:2)  = ',  mLayerTempTrial(1:2)    ! trial value of ground temperature (K)
   write(*,'(a,1x,10(f30.20))') 'scalarCanairNetNrgFlux = ', scalarCanairNetNrgFlux
   write(*,'(a,1x,10(f30.20))') 'scalarCanopyNetNrgFlux = ', scalarCanopyNetNrgFlux
   write(*,'(a,1x,10(f30.20))') 'scalarGroundNetNrgFlux = ', scalarGroundNetNrgFlux
   write(*,'(a,1x,10(f30.20))') 'dGroundNetFlux_dGroundTemp = ', dGroundNetFlux_dGroundTemp
  endif ! if checking fluxes

 endif ! if calculating the energy fluxes over vegetation

 ! *****
 ! * CALCULATE ENERGY FLUXES THROUGH THE SNOW-SOIL DOMAIN...
 ! **********************************************************

 ! check the need to compute energy fluxes throughout the snow+soil domain
 if(nSnowSoilNrg>0 .or. firstFluxCall)then

  ! calculate energy fluxes at layer interfaces through the snow and soil domain
  call ssdNrgFlux(&
                  ! input: fluxes and derivatives at the upper boundary
                  scalarGroundNetNrgFlux,                 & ! intent(in): total flux at the ground surface (W m-2)
                  dGroundNetFlux_dGroundTemp,             & ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                  ! input: liquid water fluxes throughout the snow and soil domains
                  iLayerLiqFluxSnow,                      & ! intent(in): liquid flux at the interface of each snow layer (m s-1)
                  iLayerLiqFluxSoil,                      & ! intent(in): liquid flux at the interface of each soil layer (m s-1)
                  ! input: trial value of model state variabes
                  mLayerTempTrial,                        & ! intent(in): trial temperature at the current iteration (K)
                  ! input-output: data structures
                  mpar_data,                              & ! intent(in):    model parameters
                  indx_data,                              & ! intent(in):    model indices
                  prog_data,                              & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                              & ! intent(in):    model diagnostic variables for a local HRU
                  flux_data,                              & ! intent(inout): model fluxes for a local HRU
                  ! output: fluxes and derivatives at all layer interfaces
                  iLayerNrgFlux,                          & ! intent(out): energy flux at the layer interfaces (W m-2)
                  dNrgFlux_dTempAbove,                    & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                  dNrgFlux_dTempBelow,                    & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  ! calculate net energy fluxes for each snow and soil layer (J m-3 s-1)
  do iLayer=1,nLayers
   mLayerNrgFlux(iLayer) = -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))/mLayerDepth(iLayer)
   if(globalPrintFlag)then
    if(iLayer < 3) write(*,'(a,1x,i4,1x,10(f25.15,1x))') 'iLayer, iLayerNrgFlux(iLayer-1:iLayer), mLayerNrgFlux(iLayer)   = ', iLayer, iLayerNrgFlux(iLayer-1:iLayer), mLayerNrgFlux(iLayer)
   endif
  end do

 endif  ! if computing energy fluxes throughout the snow+soil domain

 ! *****
 ! * CALCULATE THE LIQUID FLUX THROUGH VEGETATION...
 ! **************************************************

 ! check the need to compute the liquid water fluxes through vegetation
 if(ixVegHyd/=integerMissing .or. firstFluxCall)then

  ! calculate liquid water fluxes through vegetation
  call vegLiqFlux(&
                  ! input
                  computeVegFlux,                         & ! intent(in): flag to denote if computing energy flux over vegetation
                  scalarCanopyLiqTrial,                   & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                  scalarRainfall,                         & ! intent(in): rainfall rate (kg m-2 s-1)
                  ! input-output: data structures
                  mpar_data,                              & ! intent(in): model parameters
                  diag_data,                              & ! intent(in): local HRU diagnostic model variables
                  ! output
                  scalarThroughfallRain,                  & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                  scalarCanopyLiqDrainage,                & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                  scalarThroughfallRainDeriv,             & ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
                  scalarCanopyLiqDrainageDeriv,           & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  ! calculate the net liquid water flux for the vegetation canopy
  scalarCanopyNetLiqFlux = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
  
  ! calculate the total derivative in the downward liquid flux
  scalarCanopyLiqDeriv   = scalarThroughfallRainDeriv + scalarCanopyLiqDrainageDeriv
  
  ! test
  if(globalPrintFlag)then
   print*, '**'
   print*, 'scalarRainfall          = ', scalarRainfall
   print*, 'scalarThroughfallRain   = ', scalarThroughfallRain
   print*, 'scalarCanopyEvaporation = ', scalarCanopyEvaporation
   print*, 'scalarCanopyLiqDrainage = ', scalarCanopyLiqDrainage
   print*, 'scalarCanopyNetLiqFlux  = ', scalarCanopyNetLiqFlux
   print*, 'scalarCanopyLiqTrial    = ', scalarCanopyLiqTrial
  endif

 endif  ! computing the liquid water fluxes through vegetation

 ! *****
 ! * CALCULATE THE LIQUID FLUX THROUGH SNOW...
 ! ********************************************

 ! check the need to compute liquid water fluxes through snow
 if(nSnowOnlyHyd>0 .or. firstFluxCall)then

  ! compute liquid fluxes through snow
  call snowLiqFlx(&
                  ! input: model control
                  nSnow,                                 & ! intent(in): number of snow layers
                  firstFluxCall,                         & ! intent(in): the first flux call (compute variables that are constant over the iterations)
                  ! input: forcing for the snow domain
                  scalarThroughfallRain,                 & ! intent(in): rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
                  scalarCanopyLiqDrainage,               & ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
                  ! input: model state vector
                  mLayerVolFracLiqTrial(1:nSnow),        & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
                  ! input-output: data structures
                  mpar_data,                             & ! intent(in):    model parameters
                  prog_data,                             & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                             & ! intent(inout): model diagnostic variables for a local HRU
                  ! output: fluxes and derivatives
                  iLayerLiqFluxSnow(0:nSnow),            & ! intent(out): vertical liquid water flux at layer interfaces (m s-1)
                  iLayerLiqFluxSnowDeriv(0:nSnow),       & ! intent(out): derivative in vertical liquid water flux at layer interfaces (m s-1)
                  ! output: error control
                  err,cmessage)                            ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! define forcing for the soil domain
  scalarRainPlusMelt = iLayerLiqFluxSnow(nSnow)    ! drainage from the base of the snowpack

  ! calculate net liquid water fluxes for each soil layer (s-1)
  do iLayer=1,nSnow
   mLayerLiqFluxSnow(iLayer) = -(iLayerLiqFluxSnow(iLayer) - iLayerLiqFluxSnow(iLayer-1))/mLayerDepth(iLayer)
   !write(*,'(a,1x,i4,1x,2(f16.10,1x))')  'iLayer, mLayerLiqFluxSnow(iLayer), iLayerLiqFluxSnow(iLayer-1) = ', &
   !                                       iLayer, mLayerLiqFluxSnow(iLayer), iLayerLiqFluxSnow(iLayer-1)
  end do

  ! compute drainage from the soil zone (needed for mass balance checks)
  scalarSnowDrainage = iLayerLiqFluxSnow(nSnow)

 else

  ! define forcing for the soil domain for the case of no snow layers
  ! NOTE: in case where nSnowOnlyHyd==0 AND snow layers exist, then scalarRainPlusMelt is taken from the previous flux evaluation
  if(nSnow==0)then
   scalarRainPlusMelt = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water &  ! liquid flux from the canopy (m s-1)
                         + drainageMeltPond/iden_water  ! melt of the snow without a layer (m s-1)
  endif  ! if no snow layers

 endif

 ! *****
 ! * CALCULATE THE LIQUID FLUX THROUGH SOIL...
 ! ********************************************

 ! check the need to calculate the liquid flux through soil
 if(nSoilOnlyHyd>0 .or. firstFluxCall)then

  ! calculate the liquid flux through soil
  call soilLiqFlx(&
                  ! input: model control
                  nSoil,                                  & ! intent(in):    number of soil layers
                  firstSplitOper,                         & ! intent(in):    flag indicating first flux call in a splitting operation
                  .true.,                                 & ! intent(in):    flag indicating if derivatives are desired
                  ! input: trial state variables
                  mLayerTempTrial(nSnow+1:nLayers),       & ! intent(in):    trial temperature at the current iteration (K)
                  mLayerMatricHeadLiqTrial(1:nSoil),      & ! intent(in):    liquid water matric potential (m)
                  mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in):    volumetric fraction of liquid water (-)
                  mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in):    volumetric fraction of ice (-)
                  ! input: pre-computed deriavatives
                  mLayerdTheta_dTk(nSnow+1:nLayers),      & ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
                  dPsiLiq_dTemp(1:nSoil),                 & ! intent(in):    derivative in liquid water matric potential w.r.t. temperature (m K-1)
                  ! input: fluxes
                  scalarCanopyTranspiration,              & ! intent(in):    canopy transpiration (kg m-2 s-1)
                  scalarGroundEvaporation,                & ! intent(in):    ground evaporation (kg m-2 s-1)
                  scalarRainPlusMelt,                     & ! intent(in):    rain plus melt (m s-1)
                  ! input-output: data structures
                  mpar_data,                              & ! intent(in):    model parameters
                  indx_data,                              & ! intent(in):    model indices
                  prog_data,                              & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                              & ! intent(in):    model diagnostic variables for a local HRU
                  flux_data,                              & ! intent(in):    model fluxes for a local HRU
                  ! output: diagnostic variables for surface runoff
                  scalarMaxInfilRate,                     & ! intent(inout): maximum infiltration rate (m s-1)
                  scalarInfilArea,                        & ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
                  scalarFrozenArea,                       & ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
                  scalarSurfaceRunoff,                    & ! intent(out):   surface runoff (m s-1)
                  ! output: diagnostic variables for model layers
                  mLayerdTheta_dPsi,                      & ! intent(out):   derivative in the soil water characteristic w.r.t. psi (m-1)
                  mLayerdPsi_dTheta,                      & ! intent(out):   derivative in the soil water characteristic w.r.t. theta (m)
                  dHydCond_dMatric,                       & ! intent(out):   derivative in hydraulic conductivity w.r.t matric head (s-1)
                  ! output: fluxes
                  scalarInfiltration,                     & ! intent(out):   surface infiltration rate (m s-1) -- controls on infiltration only computed for iter==1
                  iLayerLiqFluxSoil,                      & ! intent(out):   liquid fluxes at layer interfaces (m s-1)
                  mLayerTranspire,                        & ! intent(out):   transpiration loss from each soil layer (m s-1)
                  mLayerHydCond,                          & ! intent(out):   hydraulic conductivity in each layer (m s-1)
                  ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                  dq_dHydStateAbove,                      & ! intent(out):   derivatives in the flux w.r.t. matric head in the layer above (s-1)
                  dq_dHydStateBelow,                      & ! intent(out):   derivatives in the flux w.r.t. matric head in the layer below (s-1)
                  ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
                  dq_dNrgStateAbove,                      & ! intent(out):   derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
                  dq_dNrgStateBelow,                      & ! intent(out):   derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 
  ! calculate net liquid water fluxes for each soil layer (s-1)
  do iLayer=1,nSoil
   mLayerLiqFluxSoil(iLayer) = -(iLayerLiqFluxSoil(iLayer) - iLayerLiqFluxSoil(iLayer-1))/mLayerDepth(iLayer+nSnow)
   !print*, 'iLayerLiqFluxSoil(iLayer-1), mLayerLiqFluxSoil(iLayer) = ', iLayerLiqFluxSoil(iLayer-1), mLayerLiqFluxSoil(iLayer)
  end do
 
  ! calculate the soil control on infiltration
  if(nSnow==0) then
   ! * case of infiltration into soil
   if(scalarMaxInfilRate > scalarRainPlusMelt)then  ! infiltration is not rate-limited
    scalarSoilControl = (1._dp - scalarFrozenArea)*scalarInfilArea
   else
    scalarSoilControl = 0._dp  ! (scalarRainPlusMelt exceeds maximum infiltration rate
   endif
  else
   ! * case of infiltration into snow
   scalarSoilControl = 1._dp
  endif
  
  ! compute drainage from the soil zone (needed for mass balance checks)
  scalarSoilDrainage = iLayerLiqFluxSoil(nSoil)

  ! expand derivatives to the total water matric potential
  ! NOTE: arrays are offset because computing derivatives in interface fluxes, at the top and bottom of the layer respectively
  if(globalPrintFlag) print*, 'dPsiLiq_dPsi0(1:nSoil) = ', dPsiLiq_dPsi0(1:nSoil)
  dq_dHydStateAbove(1:nSoil)   = dq_dHydStateAbove(1:nSoil)  *dPsiLiq_dPsi0(1:nSoil)
  dq_dHydStateBelow(0:nSoil-1) = dq_dHydStateBelow(0:nSoil-1)*dPsiLiq_dPsi0(1:nSoil)

 endif  ! if calculating the liquid flux through soil

 ! *****
 ! * CALCULATE THE GROUNDWATER FLOW...
 ! ************************************

 ! check if computing soil hydrology
 if(nSoilOnlyHyd>0 .or. firstFluxCall)then

  ! set baseflow fluxes to zero if the baseflow routine is not used
  if(local_ixGroundwater/=qbaseTopmodel)then
   ! (diagnostic variables in the data structures)
   scalarExfiltration     = 0._dp  ! exfiltration from the soil profile (m s-1)
   mLayerColumnOutflow(:) = 0._dp  ! column outflow from each soil layer (m3 s-1)
   ! (variables needed for the numerical solution)
   mLayerBaseflow(:)      = 0._dp  ! baseflow from each soil layer (m s-1)
  
  ! topmodel-ish shallow groundwater
  else ! local_ixGroundwater==qbaseTopmodel
  
   ! check the derivative matrix is sized appropriately
   if(size(dBaseflow_dMatric,1)/=nSoil .or. size(dBaseflow_dMatric,2)/=nSoil)then
    message=trim(message)//'expect dBaseflow_dMatric to be nSoil x nSoil'
    err=20; return
   endif
 
   ! compute the baseflow flux
   call groundwatr(&
                   ! input: model control
                   nSnow,                                   & ! intent(in):    number of snow layers
                   nSoil,                                   & ! intent(in):    number of soil layers
                   nLayers,                                 & ! intent(in):    total number of layers
                   firstFluxCall,                           & ! intent(in):    logical flag to compute index of the lowest saturated layer
                   ! input: state and diagnostic variables
                   mLayerdTheta_dPsi,                       & ! intent(in):    derivative in the soil water characteristic w.r.t. matric head in each layer (m-1)
                   mLayerMatricHeadLiqTrial,                & ! intent(in):    liquid water matric potential (m)
                   mLayerVolFracLiqTrial(nSnow+1:nLayers),  & ! intent(in):    volumetric fraction of liquid water (-)
                   mLayerVolFracIceTrial(nSnow+1:nLayers),  & ! intent(in):    volumetric fraction of ice (-)
                   ! input: data structures
                   attr_data,                               & ! intent(in):    model attributes
                   mpar_data,                               & ! intent(in):    model parameters
                   prog_data,                               & ! intent(in):    model prognostic variables for a local HRU
                   diag_data,                               & ! intent(in):    model diagnostic variables for a local HRU
                   flux_data,                               & ! intent(inout): model fluxes for a local HRU
                   ! output
                   ixSaturation,                            & ! intent(inout) index of lowest saturated layer (NOTE: only computed on the first iteration)
                   mLayerBaseflow,                          & ! intent(out): baseflow from each soil layer (m s-1)
                   dBaseflow_dMatric,                       & ! intent(out): derivative in baseflow w.r.t. matric head (s-1)
                   err,cmessage)                              ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif  ! computing baseflow flux
  
  ! compute total baseflow from the soil zone (needed for mass balance checks)
  scalarSoilBaseflow = sum(mLayerBaseflow)

 endif  ! if computing soil hydrology

 ! *****
 ! (7) CALCUALTE FLUXES FOR THE DEEP AQUIFER...
 ! ********************************************

 ! identify modeling decision
 if(local_ixGroundwater==bigBucket)then
  ! deep aquifer is not yet transferred from old code structure
  message=trim(message)//'bigBucket groundwater parameterization is not yet transfered from old code structure'
  err=20; return
 else
  ! if no aquifer, then fluxes are zero
  scalarAquiferTranspire = 0._dp  ! transpiration loss from the aquifer (m s-1
  scalarAquiferRecharge  = 0._dp  ! recharge to the aquifer (m s-1)
  scalarAquiferBaseflow  = 0._dp  ! total baseflow from the aquifer (m s-1)
 end if

 ! *****
 ! (X) WRAP UP...
 ! *************
 ! define model flux vector for the vegetation sub-domain
 if(ixCasNrg/=integerMissing) fluxVec(ixCasNrg) = scalarCanairNetNrgFlux/canopyDepth
 if(ixVegNrg/=integerMissing) fluxVec(ixVegNrg) = scalarCanopyNetNrgFlux/canopyDepth
 if(ixVegHyd/=integerMissing) fluxVec(ixVegHyd) = scalarCanopyNetLiqFlux   ! NOTE: solid fluxes are handled separately

 ! populate the flux vector for energy 
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   fluxVec( ixSnowSoilNrg(iLayer) ) = mLayerNrgFlux(iLayer)
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! populate the flux vector for hydrology 
 ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
 if(nSnowSoilHyd>0)then  ! check if any hydrology states exist
  do iLayer=1,nLayers
   if(ixSnowSoilHyd(iLayer)/=integerMissing)then   ! check if a given hydrology state exists 
    select case( layerType(iLayer) )
     case(iname_snow); fluxVec( ixSnowSoilHyd(iLayer) ) = mLayerLiqFluxSnow(iLayer)
     case(iname_soil); fluxVec( ixSnowSoilHyd(iLayer) ) = mLayerLiqFluxSoil(iLayer-nSnow)
     case default; err=20; message=trim(message)//'expect layerType to be either iname_snow or iname_soil'; return
    end select
   endif  ! if a given hydrology state exists
  end do ! looping through non-missing energy state variables in the snow+soil domain
 endif  ! if any hydrology states exist

 ! set the first flux call to false
 firstFluxCall=.false.

 ! end association to variables in the data structures
 end associate

 end subroutine computFlux


 ! **********************************************************************************************************
 ! public subroutine soilCmpres: compute soil compressibility (-) and its derivative w.r.t matric head (m-1)
 ! **********************************************************************************************************
 subroutine soilCmpres(&
                       ! input:
                       ixRichards,                         & ! intent(in): choice of option for Richards' equation
                       mLayerMatricHead,                   & ! intent(in): matric head at the start of the time step (m)
                       mLayerMatricHeadTrial,              & ! intent(in): trial value of matric head (m)
                       mLayerVolFracLiqTrial,              & ! intent(in): trial value for the volumetric liquid water content in each soil layer (-)
                       mLayerVolFracIceTrial,              & ! intent(in): trial value for the volumetric ice content in each soil layer (-)
                       mLayerdTheta_dPsi,                  & ! intent(in): derivative in the soil water characteristic (m-1)
                       specificStorage,                    & ! intent(in): specific storage coefficient (m-1)
                       theta_sat,                          & ! intent(in): soil porosity (-)
                       ! output:
                       compress,                           & ! intent(out): compressibility of the soil matrix (-)
                       dCompress_dPsi,                     & ! intent(out): derivative in compressibility w.r.t. matric head (m-1)
                       err,message)                          ! intent(out): error code and error message
 implicit none
 ! input:
 integer(i4b),intent(in)        :: ixRichards                ! choice of option for Richards' equation
 real(dp),intent(in)            :: mLayerMatricHead(:)       ! matric head at the start of the time step (m)
 real(dp),intent(in)            :: mLayerMatricHeadTrial(:)  ! trial value for matric head (m)
 real(dp),intent(in)            :: mLayerVolFracLiqTrial(:)  ! trial value for volumetric fraction of liquid water (-)
 real(dp),intent(in)            :: mLayerVolFracIceTrial(:)  ! trial value for volumetric fraction of ice (-)
 real(dp),intent(in)            :: mLayerdTheta_dPsi(:)      ! derivative in the soil water characteristic (m-1)
 real(dp),intent(in)            :: specificStorage           ! specific storage coefficient (m-1)
 real(dp),intent(in)            :: theta_sat(:)              ! soil porosity (-)
 ! output:
 real(dp),intent(out)           :: compress(:)               ! soil compressibility (-)
 real(dp),intent(out)           :: dCompress_dPsi(:)         ! derivative in soil compressibility w.r.t. matric head (m-1)
 integer(i4b),intent(out)       :: err                       ! error code
 character(*),intent(out)       :: message                   ! error message
 ! local variables
 real(dp)                       :: volFracWat                ! total volumetric fraction of water (-)
 real(dp)                       :: fPart1,fPart2             ! different parts of the function
 real(dp)                       :: dPart1,dPart2             ! derivatives for different parts of the function
 integer(i4b)                   :: iLayer                    ! index of soil layer
 ! --------------------------------------------------------------
 ! initialize error control
 err=0; message='soilCmpres/'
 ! (only compute for the mixed form of Richards' equation)
 if(ixRichards==mixdform)then
  do iLayer=1,size(mLayerMatricHead)
   ! compute the total volumetric fraction of water (-)
   volFracWat = mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer)
   ! compute the compressibility term (-)
   compress(iLayer) = (specificStorage*volFracWat/theta_sat(iLayer)) * (mLayerMatricHeadTrial(iLayer) - mLayerMatricHead(iLayer))
   ! compute the derivative for the compressibility term (m-1)
   fPart1 = specificStorage*(volFracWat/theta_sat(iLayer))  ! function for the 1st part (m-1)
   fPart2 = mLayerMatricHeadTrial(iLayer) - mLayerMatricHead(iLayer)     ! function for the 2nd part (m)
   dPart1 = mLayerdTheta_dPsi(iLayer)*specificStorage/theta_sat(iLayer)  ! derivative for the 1st part (m-2)
   dPart2 = 1._dp                                                        ! derivative for the 2nd part (-)
   dCompress_dPsi(iLayer) = fPart1*dPart2 + dPart1*fPart2                ! m-1
  end do
 else
  compress(:)       = 0._dp
  dCompress_dPsi(:) = 0._dp
 end if
 end subroutine soilCmpres

end module computFlux_module
