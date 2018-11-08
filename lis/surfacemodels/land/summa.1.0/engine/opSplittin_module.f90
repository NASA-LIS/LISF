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

module opSplittin_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access matrix information
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

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

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookATTR       ! named variables for structure elements
USE var_lookup,only:iLookTYPE       ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

! provide access to the number of flux variables
USE var_lookup,only:nFlux=>maxvarFlux ! number of model flux variables

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
USE mDecisions_module,only:      &
 qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                      & ! a big bucket (lumped aquifer model)
 noExplicit                        ! no explicit groundwater parameterization

! safety: set private unless specified otherwise
implicit none
private
public::opSplittin

! named variables for the solution method
integer(i4b),parameter  :: fullyCoupled=1             ! 1st try: fully coupled solution
integer(i4b),parameter  :: splitStateType=2           ! 2nd try: split the solution by state type (energy and water)
integer(i4b),parameter  :: splitDomainType=3          ! 3rd try: split the solution by domain type (veg, snow, and soil)
integer(i4b),parameter  :: explicitEuler=4            ! 4th try: explicit Euler solution for sub-domains of a given type

! named variables for the state variable split
integer(i4b),parameter  :: nrgSplit=1                 ! order in sequence for the energy operation
integer(i4b),parameter  :: massSplit=2                ! order in sequence for the mass operation

! named variables for the domain type split
integer(i4b),parameter  :: vegSplit=1                 ! order in sequence for the vegetation split
integer(i4b),parameter  :: snowSplit=2                ! order in sequence for the snow split
integer(i4b),parameter  :: soilSplit=3                ! order in sequence for the soil split

! maximum number of possible splits
integer(i4b),parameter  :: nStateTypes=2              ! number of state types (energy, water)
integer(i4b),parameter  :: nDomains=3                 ! number of domains (vegetation, snow, and soil)

! control parameters
real(dp),parameter      :: valueMissing=-9999._dp     ! missing value
real(dp),parameter      :: verySmall=1.e-12_dp        ! a very small number (used to check consistency)
real(dp),parameter      :: veryBig=1.e+20_dp          ! a very big number
real(dp),parameter      :: dx = 1.e-8_dp              ! finite difference increment

contains


 ! **********************************************************************************************************
 ! public subroutine opSplittin: run the coupled energy-mass model for one timestep
 !
 ! The logic of the solver is as follows:
 ! (1) Attempt different solutions in the following order: (a) fully coupled; (b) split by state type (energy
 !      and mass); (c) split by domain type or a given energy and mass split (vegetation, snow, and soil);
 !      and (d) explicit Euler solution for a given state type and domain subset.
 ! (2) For a given split, compute a variable number of substeps (in varSubstep).
 ! **********************************************************************************************************
 subroutine opSplittin(&
                       ! input: model control
                       nSnow,          & ! intent(in):    number of snow layers
                       nSoil,          & ! intent(in):    number of soil layers
                       nLayers,        & ! intent(in):    total number of layers
                       nState,         & ! intent(in):    total number of state variables
                       dt,             & ! intent(inout): time step (s)
                       firstSubStep,   & ! intent(in):    flag to denote first sub-step
                       computeVegFlux, & ! intent(in):    flag to denote if computing energy flux over vegetation
                       ! input/output: data structures
                       type_data,      & ! intent(in):    type of vegetation and soil
                       attr_data,      & ! intent(in):    spatial attributes
                       forc_data,      & ! intent(in):    model forcing data
                       mpar_data,      & ! intent(in):    model parameters
                       indx_data,      & ! intent(inout): index data
                       prog_data,      & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,      & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,      & ! intent(inout): model fluxes for a local HRU
                       bvar_data,      & ! intent(in):    model variables for the local basin
                       model_decisions,& ! intent(in):    model decisions
                       ! output: model control
                       dtMultiplier,   & ! intent(out):   substep multiplier (-)
                       tooMuchMelt,    & ! intent(out):   flag to denote that ice is insufficient to support melt
                       stepFailure,    & ! intent(out):   flag to denote step failure
                       ixSolution,     & ! intent(out):   solution method used in this iteration
                       err,message)      ! intent(out):   error code and error message
 ! ---------------------------------------------------------------------------------------
 ! structure allocations
 USE globalData,only:flux_meta                        ! metadata on the model fluxes
 USE globalData,only:diag_meta                        ! metadata on the model diagnostic variables
 USE globalData,only:prog_meta                        ! metadata on the model prognostic variables
 USE globalData,only:deriv_meta                       ! metadata on the model derivatives
 USE globalData,only:flux2state_orig                  ! metadata on flux-to-state mapping (original state variables)
 USE globalData,only:flux2state_liq                   ! metadata on flux-to-state mapping (liquid water state variables)
 USE allocspace_module,only:allocLocal                ! allocate local data structures
 ! simulation of fluxes and residuals given a trial state vector
 USE soil_utils_module,only:matricHead                ! compute the matric head based on volumetric water content
 USE soil_utils_module,only:liquidHead                ! compute the liquid water matric potential
 ! population/extraction of state vectors
 USE indexState_module,only:indexSplit                ! get state indices
 USE varSubstep_module,only:varSubstep                ! complete substeps for a given split
 ! numerical recipes utility modules
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 integer(i4b),intent(in)         :: nSnow                          ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                          ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                        ! total number of layers
 integer(i4b),intent(in)         :: nState                         ! total number of state variables
 real(dp),intent(inout)          :: dt                             ! time step (seconds)
 logical(lgt),intent(in)         :: firstSubStep                   ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)         :: computeVegFlux                 ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input/output: data structures
 type(var_i),intent(in)          :: type_data                      ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                      ! spatial attributes
 type(var_d),intent(in)          :: forc_data                      ! model forcing data
 type(var_dlength),intent(in)    :: mpar_data                      ! model parameters
 type(var_ilength),intent(inout) :: indx_data                      ! indices for a local HRU
 type(var_dlength),intent(inout) :: prog_data                      ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                      ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                      ! model fluxes for a local HRU
 type(var_dlength),intent(in)    :: bvar_data                      ! model variables for the local basin
 type(model_options),intent(in)  :: model_decisions(:)             ! model decisions
 ! output: model control
 real(dp),intent(out)            :: dtMultiplier                   ! substep multiplier (-)
 logical(lgt),intent(out)        :: tooMuchMelt                    ! flag to denote that ice is insufficient to support melt
 logical(lgt),intent(out)        :: stepFailure                    ! flag to denote step failure 
 integer(i4b),intent(out)        :: err                            ! error code
 character(*),intent(out)        :: message                        ! error message
 ! *********************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************
 ! ---------------------------------------------------------------------------------------
 ! * general local variables
 ! ---------------------------------------------------------------------------------------
 character(LEN=256)              :: cmessage                       ! error message of downwind routine
 integer(i4b)                    :: iSoil                          ! index of soil layer
 integer(i4b)                    :: iVar                           ! index of variables in data structures
 logical(lgt)                    :: firstSuccess                   ! flag to define the first success
 logical(lgt)                    :: firstFluxCall                  ! flag to define the first flux call
 logical(lgt)                    :: reduceCoupledStep              ! flag to define the need to reduce the length of the coupled step
 type(var_dlength)               :: prog_temp                      ! temporary model prognostic variables
 type(var_dlength)               :: diag_temp                      ! temporary model diagnostic variables
 type(var_dlength)               :: deriv_data                     ! derivatives in model fluxes w.r.t. relevant state variables 
 real(dp),dimension(nLayers)     :: mLayerVolFracIceInit           ! initial vector for volumetric fraction of ice (-)
 ! ------------------------------------------------------------------------------------------------------
 ! * operator splitting
 ! ------------------------------------------------------------------------------------------------------
 ! minimum time step
 real(dp)                        :: dt_min                         ! minimum time step (seconds)
 real(dp),parameter              :: dtmin_fullyCoupled=10._dp      ! minimum time step for the fully coupled solution
 real(dp),parameter              :: dtmin_splitStateType=1._dp     ! minimum time step for the split by state type 
 real(dp),parameter              :: dtmin_splitDomainType=0.1_dp   ! minimum time step for the split by domain type 
 real(dp),parameter              :: dtmin_explicitEuler=0.1_dp     ! minimum time step for the explicit Euler solution
 ! explicit error tolerance (depends on state type split, so defined here)
 real(dp),parameter              :: errorTolLiqFlux=0.01_dp        ! error tolerance in the explicit solution (liquid flux)
 real(dp),parameter              :: errorTolNrgFlux=10._dp         ! error tolerance in the explicit solution (energy flux)
 real(dp)                        :: errTol                         ! error tolerance in the explicit solution
 ! number of substeps taken for a given split
 integer(i4b)                    :: nSubsteps                      ! number of substeps taken for a given split
 ! named variables defining the solution method
 integer(i4b)                    :: ixSolution                     ! index of solution method (1,2,3,...)
 ! actual number of splits
 integer(i4b)                    :: nStateTypeSplit                ! number of splits for the state type
 integer(i4b)                    :: nDomainSplit                   ! number of splits for the domain
 ! indices for the state type split
 integer(i4b)                    :: iStateTypeSplit                ! index of the state type split
 integer(i4b)                    :: iTrialStateSplit               ! index of state split trial
 ! indices for the domain split
 integer(i4b)                    :: iDomainSplit                   ! index of the domain split
 ! state and flux masks for a given split
 integer(i4b),dimension(nState)  :: stateCheck                     ! number of times each state variable is updated (should=1) 
 logical(lgt),dimension(nState)  :: stateMask                      ! mask defining desired state variables
 logical(lgt),dimension(nFlux)   :: fluxMask                       ! mask defining desired flux variables
 integer(i4b)                    :: nSubset                        ! number of selected state variables for a given split
 ! flags
 logical(lgt)                    :: failure                        ! flag to denote failure of substepping
 logical(lgt)                    :: doAdjustTemp                   ! flag to adjust temperature after the mass split
 logical(lgt)                    :: failedMinimumStep              ! flag to denote failure of substepping for a given split
 integer(i4b)                    :: ixSaturation                   ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 globalVars: associate(&
 ! model decisions
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,& ! intent(in):    [i4b]    groundwater parameterization
 ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,& ! intent(in):    [i4b]    spatial representation of groundwater (local-column or single-basin)
 ! domain boundary conditions
 airtemp                 => forc_data%var(iLookFORCE%airtemp)                      ,& ! intent(in):    [dp]     temperature of the upper boundary of the snow and soil domains (K)
 ! vector of energy and hydrology indices for the snow and soil domains
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):    [i4b]    number of hydrology state variables in the snow+soil domain
 ! indices of model state variables
 ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in):    [i4b(:)] indices defining the type of the state (ixNrgState...)
 ixNrgCanair             => indx_data%var(iLookINDEX%ixNrgCanair)%dat              ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for energy states in canopy air space domain
 ixNrgCanopy             => indx_data%var(iLookINDEX%ixNrgCanopy)%dat              ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for energy states in the canopy domain
 ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
 ixNrgLayer              => indx_data%var(iLookINDEX%ixNrgLayer)%dat               ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer              => indx_data%var(iLookINDEX%ixHydLayer)%dat               ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy energy state variable
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)
 ! domain configuration
 canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in):    [dp]     canopy depth (m)
 mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):    [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
 ! snow parameters
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)         ,& ! intent(in):    [dp]     scaling parameter for the snow freezing curve (K-1)
 ! depth-varying soil parameters
 vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat               ,& ! intent(in):    [dp(:)]  van Genutchen "m" parameter (-)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)%dat                    ,& ! intent(in):    [dp(:)]  van Genutchen "n" parameter (-)
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)%dat                ,& ! intent(in):    [dp(:)]  van Genutchen "alpha" parameter (m-1)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,& ! intent(in):    [dp(:)]  soil porosity (-)
 theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat                ,& ! intent(in):    [dp(:)]  soil residual volumetric water content (-)
 ! soil parameters
 specificStorage         => mpar_data%var(iLookPARAM%specificStorage)%dat(1)       ,& ! intent(in):    [dp]     specific storage coefficient (m-1)
 ! model diagnostic variables (fraction of liquid water)
 scalarFracLiqVeg        => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)       ,& ! intent(out):   [dp]     fraction of liquid water on vegetation (-)
 mLayerFracLiqSnow       => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat         ,& ! intent(out):   [dp(:)]  fraction of liquid water in each snow layer (-)
 mLayerMeltFreeze        => diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat          ,& ! intent(out):   [dp(:)]  melt-freeze in each snow and soil layer (kg m-3)
 ! model state variables (vegetation canopy)
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(out):   [dp]     temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(out):   [dp]     temperature of the vegetation canopy (K)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(out):   [dp]     mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(out):   [dp]     mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(out):   [dp]     mass of total water on the vegetation canopy (kg m-2)
 ! model state variables (snow and soil domains)
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(out):   [dp(:)]  temperature of each snow/soil layer (K)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(out):   [dp(:)]  volumetric fraction of ice (-)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(out):   [dp(:)]  volumetric fraction of liquid water (-)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(out):   [dp(:)]  volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(out):   [dp(:)]  matric head (m)
 mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat        & ! intent(out):   [dp(:)]  matric potential of liquid water (m)
 )
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="opSplittin/"

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! -----
 ! * initialize...
 ! ---------------

 ! set the global print flag
 globalPrintFlag=.false.

 if(globalPrintFlag)&
 print*, trim(message), dt

 ! initialize the first flux call
 firstSuccess=.false.
 firstFluxCall=.true.

 ! initialize the flags 
 tooMuchMelt=.false.  ! too much melt (merge snow layers)
 stepFailure=.false.  ! step failure

 ! initialize flag for the success of the substepping
 failure=.false.

 ! initialize the state check
 stateCheck(:) = 0

 ! compute the total water content in the vegetation canopy
 scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce  ! kg m-2

 ! save volumetric ice content at the start of the step
 ! NOTE: used for volumetric loss due to melt-freeze
 mLayerVolFracIceInit(:) = mLayerVolFracIce(:)

 ! compute the total water content in snow and soil
 ! NOTE: no ice expansion allowed for soil
 if(nSnow>0)& 
 mLayerVolFracWat(      1:nSnow  ) = mLayerVolFracLiq(      1:nSnow  ) + mLayerVolFracIce(      1:nSnow  )*(iden_ice/iden_water)
 mLayerVolFracWat(nSnow+1:nLayers) = mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)

 ! compute the liquid water matric potential (m)
 ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
 ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
 do iSoil=1,nSoil
  call liquidHead(mLayerMatricHead(iSoil),mLayerVolFracLiq(nSnow+iSoil),mLayerVolFracIce(nSnow+iSoil),  & ! input:  state variables
                  vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m(iSoil),         & ! input:  parameters
                  matricHeadLiq=mLayerMatricHeadLiq(iSoil),                                             & ! output: liquid water matric potential (m)
                  err=err,message=cmessage)                                                               ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through soil layers (computing liquid water matric potential)

 ! allocate space for the temporary prognostic variable structure
 call allocLocal(prog_meta(:),prog_temp,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the temporary diagnostic variable structure
 call allocLocal(diag_meta(:),diag_temp,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the derivative structure
 call allocLocal(deriv_meta(:),deriv_data,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! initialize solution method
 ixSolution=fullyCoupled

 ! loop through solution methods 
 solution: do
 
  ! define the number of operator splits for the state type
  if(ixSolution==fullyCoupled)then
   nStateTypeSplit=1
  else
   nStateTypeSplit=nStateTypes
  endif

  ! state type splitting loop
  stateTypeSplit: do iStateTypeSplit=1,nStateTypeSplit 
   trialStateSplit: do iTrialStateSplit=1,2  ! two state type trials: 1=full domain; 2=sub-domains
  
    ! define the number of operator splits for the domain
    if(ixSolution==fullyCoupled .or. ixSolution==splitStateType)then
     nDomainSplit=1
    else
     nDomainSplit=nDomains
    endif
  
    ! flag to adjust the temperature
    doAdjustTemp = (ixSolution/=fullyCoupled .and. iStateTypeSplit==massSplit)
  
    ! get the error tolerance
    errTol = merge(errorTolNrgFlux,errorTolLiqFlux,iStateTypeSplit==nrgSplit)
  
    ! -----
    ! * modify state variables for the mass split...
    ! ----------------------------------------------
   
    ! modify the state type names associated with the state vector
    if(ixSolution/=fullyCoupled .and. iStateTypeSplit==massSplit)then
     if(computeVegFlux)then
      where(ixStateType(ixHydCanopy)==iname_watCanopy) ixStateType(ixHydCanopy)=iname_liqCanopy
     endif
     where(ixStateType(ixHydLayer) ==iname_watLayer)  ixStateType(ixHydLayer) =iname_liqLayer
     where(ixStateType(ixHydLayer) ==iname_matLayer)  ixStateType(ixHydLayer) =iname_lmpLayer
    endif  ! if modifying state variables for the mass split
   
    ! domain type splitting loop
    domainSplit: do iDomainSplit=1,nDomainSplit
  
     ! try multiple solution methods
     trySolution: do    ! exit upon success
   
      ! *******************************************************************************************************************************
      ! *******************************************************************************************************************************
      ! *******************************************************************************************************************************
      ! ***** trial with a given solution method...
     
      ! initialize error control
      err=0; message="opSplittin/"

      ! -----
      ! * define subsets for a given split...
      ! -------------------------------------

      ! get the mask for the state subset
      call stateFilter(ixSolution,iStateTypeSplit,iDomainSplit,indx_data,stateMask,nSubset,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

      ! check
      if(globalPrintFlag)then
       print*, 'after filter: stateMask   = ', stateMask
       print*, 'ixSolution==explicitEuler = ', ixSolution==explicitEuler
      endif

      ! check that state variables exist
      if(nSubset==0) cycle domainSplit 
 
      ! -----
      ! * assemble vectors for a given split...
      ! ---------------------------------------
    
      ! define minimum time step
      select case(ixSolution)
       case(fullyCoupled);    dt_min = dtmin_fullyCoupled
       case(splitStateType);  dt_min = dtmin_splitStateType
       case(splitDomainType); dt_min = dtmin_splitDomainType
       case(explicitEuler);   dt_min = dtmin_explicitEuler
       case default; err=20; message=trim(message)//'solution method not found'; return
      end select

      ! get indices for a given split
      call indexSplit(stateMask,                   & ! intent(in)    : logical vector (.true. if state is in the subset)
                      nSnow,nSoil,nLayers,nSubset, & ! intent(in)    : number of snow and soil layers, and total number of layers
                      indx_data,                   & ! intent(inout) : index data structure
                      err,cmessage)                  ! intent(out)   : error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
     
      ! define the mask of the fluxes used
      stateSubset: associate(ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat)
      do iVar=1,size(flux_meta)

       ! * split solution
       if(ixSolution/=fullyCoupled)then
        select case(iStateTypeSplit)
         case(nrgSplit);  fluxMask(iVar) = any(ixStateType_subset==flux2state_orig(iVar)%state1) .or. any(ixStateType_subset==flux2state_orig(iVar)%state2) 
         case(massSplit); fluxMask(iVar) = any(ixStateType_subset==flux2state_liq(iVar)%state1)  .or. any(ixStateType_subset==flux2state_liq(iVar)%state2)
         case default; err=20; message=trim(message)//'unable to identify split based on state type'; return
        end select
  
       ! * fully coupled
       else
        fluxMask(iVar) = any(ixStateType_subset==flux2state_orig(iVar)%state1) .or. any(ixStateType_subset==flux2state_orig(iVar)%state2)
       endif
       
       ! * check
       if(globalPrintFlag .and. fluxMask(iVar))&
       print*, trim(flux_meta(iVar)%varname)
  
      end do
      end associate stateSubset
   
      ! initialize the model fluxes (some model fluxes are not computed in the iterations)
      do iVar=1,size(flux_meta)
       if(fluxMask(iVar)) flux_data%var(iVar)%dat(:) = 0._dp
      end do
   
      ! -----
      ! * solve variable subset for one time step...
      ! --------------------------------------------

      ! reset the flag for the first flux call
      if(.not.firstSuccess) firstFluxCall=.true.
 
      ! save/recover copies of prognostic variables
      do iVar=1,size(prog_data%var)
       select case(failure)
        case(.false.); prog_temp%var(iVar)%dat(:) = prog_data%var(iVar)%dat(:)
        case(.true.);  prog_data%var(iVar)%dat(:) = prog_temp%var(iVar)%dat(:)
       end select
      end do  ! looping through variables
   
      ! save/recover copies of diagnostic variables
      do iVar=1,size(diag_data%var)
       select case(failure)
        case(.false.); diag_temp%var(iVar)%dat(:) = diag_data%var(iVar)%dat(:)
        case(.true.);  diag_data%var(iVar)%dat(:) = diag_temp%var(iVar)%dat(:)
       end select
      end do  ! looping through variables

      ! solve variable subset for one full time step
      call varSubstep(&
                      ! input: model control
                      dt,                         & ! intent(inout) : time step (s)
                      dt_min,                     & ! intent(in)    : minimum time step (seconds)
                      errTol,                     & ! intent(in)    : error tolerance for the explicit solution
                      nSubset,                    & ! intent(in)    : total number of variables in the state subset
                      doAdjustTemp,               & ! intent(in)    : flag to indicate if we adjust the temperature
                      firstSubStep,               & ! intent(in)    : flag to denote first sub-step
                      firstFluxCall,              & ! intent(inout) : flag to indicate if we are processing the first flux call
                      (ixSolution==explicitEuler),& ! intent(in)    : flag to denote computing the explicit Euler solution
                      computeVegFlux,             & ! intent(in)    : flag to denote if computing energy flux over vegetation
                      fluxMask,                   & ! intent(in)    : mask for the fluxes used in this given state subset
                      ! input/output: data structures
                      model_decisions,            & ! intent(in)    : model decisions
                      type_data,                  & ! intent(in)    : type of vegetation and soil
                      attr_data,                  & ! intent(in)    : spatial attributes
                      forc_data,                  & ! intent(in)    : model forcing data
                      mpar_data,                  & ! intent(in)    : model parameters
                      indx_data,                  & ! intent(inout) : index data
                      prog_data,                  & ! intent(inout) : model prognostic variables for a local HRU
                      diag_data,                  & ! intent(inout) : model diagnostic variables for a local HRU
                      flux_data,                  & ! intent(inout) : model fluxes for a local HRU
                      deriv_data,                 & ! intent(inout) : derivatives in model fluxes w.r.t. relevant state variables
                      bvar_data,                  & ! intent(in)    : model variables for the local basin
                      ! output: control
                      ixSaturation,               & ! intent(inout) : index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      dtMultiplier,               & ! intent(out)   : substep multiplier (-)
                      nSubsteps,                  & ! intent(out)   : number of substeps taken for a given split
                      failedMinimumStep,          & ! intent(out)   : flag for failed substeps
                      reduceCoupledStep,          & ! intent(out)   : flag to reduce the length of the coupled step
                      tooMuchMelt,                & ! intent(out)   : flag to denote that ice is insufficient to support melt
                      err,cmessage)                 ! intent(out)   : error code and error message
      if(err/=0)then
       message=trim(message)//trim(cmessage)
       if(err>0) return
      endif  ! (check for errors)

      ! check 
      if(globalPrintFlag .and. ixSolution>splitStateType)then
       print*, 'dt = ', dt
       print*, 'after varSubstep: err            = ', err
       print*, 'after varSubstep: cmessage       = ', trim(cmessage)
       print*, 'after varSubstep: stateMask      = ', stateMask
       print*, 'iStateTypeSplit, nStateTypeSplit = ', iStateTypeSplit, nStateTypeSplit
       print*, 'iDomainSplit,    nDomainSplit    = ', iDomainSplit,    nDomainSplit
       print*, 'nSubset           = ', nSubset
       print*, 'tooMuchMelt       = ', tooMuchMelt
       print*, 'reduceCoupledStep = ', reduceCoupledStep
       print*, 'failedMinimumStep = ', failedMinimumStep, merge('coupled','opSplit',ixSolution==fullyCoupled)
       !if(ixSolution==explicitEuler)then
       ! print*, trim(message)//trim(cmessage)
       ! print*, 'PAUSE: failed splitStateType attempt'; read(*,*)
       !endif
      endif    

      ! if too much melt then return
      ! NOTE: need to go all the way back to coupled_em and merge snow layers, as all splitting operations need to occur with the same layer geometry
      if(tooMuchMelt .or. reduceCoupledStep)then
       stepFailure=.true.
       return
      endif
 
      ! define failure
      failure = (failedMinimumStep .or. err<0)
      if(.not.failure) firstSuccess=.true. 
 
      ! -----
      ! * success: revert back to "more coupled" methods...
      ! ---------------------------------------------------
     
      ! success = exit the trySolution loop
      if(.not.failure)then

       ! check that state variables updated
       where(stateMask) stateCheck = stateCheck+1
       if(any(stateCheck>1))then
        message=trim(message)//'state variable updated more than once!'
        err=20; return
       endif

       ! fully coupled
       if(ixSolution==fullyCoupled)then
        exit stateTypeSplit     ! stay within the ixSolution loop in case mass balance errors

       ! split state type
       elseif(ixSolution==splitStateType)then
        exit domainSplit  ! this (1) tries the other state split if iStateTypeSplit<nStateTypeSplit, or (2) exits since finished

       ! split domain type
       elseif(ixSolution==splitDomainType)then
        if(iDomainSplit==nDomainSplit)then
         ixSolution=splitStateType   ! try to solve all domains simultaneously for the next state split
         exit domainSplit            
        else
         exit trySolution            ! iDomainSplit<nDomainSplit: try the next domain
        endif

       ! explicit Euler
       elseif(ixSolution==explicitEuler)then
        if(iDomainSplit==nDomainSplit)then
         ixSolution=splitStateType   ! try to solve all domains simultaneously for the next state split (implicit solution)
         exit domainSplit
        else
         ixSolution=splitDomainType  ! try implicit solution for the next domain
         exit trySolution            ! iDomainSplit<nDomainSplit: try the next domain
        endif

       ! check that the index of ixSolution is known
       else
        message=trim(message)//'unknown solution index'
        err=20; return
       endif

      ! -----
      ! * if failed substeps then try another solution method...
      ! --------------------------------------------------------
     
      ! adjust solution if failed to converge at the prescribed minimum time step dt_min
      else
   
       ! if fully coupled then start again by splitting by state type
       if(ixSolution==fullyCoupled)then
        ixSolution=splitStateType
        exit stateTypeSplit  ! split into state types
   
       ! if splitting by state type then start again by splitting by domain type
       ! NOTE: will do another trial using the trialStateSplit loop 
       elseif(ixSolution==splitStateType)then
        ixSolution=splitDomainType
        exit domainSplit
   
       ! if already splitting by domain type then move to the last resort of explicit Euler
       elseif(ixSolution==splitDomainType)then
        ixSolution=explicitEuler
        cycle trySolution
   
       ! check that did not fail with explicit Euler (this was the last resort!)
       elseif(ixSolution==explicitEuler)then
        stepFailure=.true.
        return
   
       ! check that the index of ixSolution is known
       else
        message=trim(message)//'unknown solution index'
        err=20; return
       endif
        
      endif  ! if failed substeps
   
      ! ***** trial with a given solution method...
      ! *******************************************************************************************************************************
      ! *******************************************************************************************************************************
      ! *******************************************************************************************************************************
  
     end do trySolution  ! trial of given solution method
    end do domainSplit ! domain type splitting loop
  
    ! -----
    ! * reset state variables for the mass split...
    ! ---------------------------------------------
   
    ! modify the state type names associated with the state vector
    if(ixSolution/=fullyCoupled .and. iStateTypeSplit==massSplit)then
     if(computeVegFlux)then
      where(ixStateType(ixHydCanopy)==iname_liqCanopy) ixStateType(ixHydCanopy)=iname_watCanopy
     endif
     where(ixStateType(ixHydLayer) ==iname_liqLayer)  ixStateType(ixHydLayer) =iname_watLayer
     where(ixStateType(ixHydLayer) ==iname_lmpLayer)  ixStateType(ixHydLayer) =iname_matLayer
    endif  ! if modifying state variables for the mass split

    ! success = exit the trialStateSplit loop
    if(.not.failure) exit trialStateSplit

   end do trialStateSplit ! trial different state type splits  
  end do stateTypeSplit ! state type splitting loop 
  !print*, 'PAUSE: end of splitting loop'; read(*,*)

  ! ==========================================================================================================================================
  ! ==========================================================================================================================================

  ! success = exit the trySolution loop
  if(.not.failure) exit solution

 end do solution ! solution method

 ! check that all state variables were updated
 if(any(stateCheck==0))then
  message=trim(message)//'some state variables were not updated!'
  err=20; return
 endif

 ! use step halving if unable to complete the fully coupled solution in one substep
 if(ixSolution/=fullyCoupled .or. nSubsteps>1) dtMultiplier=0.5_dp

 ! compute the melt in each snow and soil layer
 if(nSnow>0) mLayerMeltFreeze(      1:nSnow  ) = -(mLayerVolFracIce(      1:nSnow  ) - mLayerVolFracIceInit(      1:nSnow  ))*iden_ice
             mLayerMeltFreeze(nSnow+1:nLayers) = -(mLayerVolFracIce(nSnow+1:nLayers) - mLayerVolFracIceInit(nSnow+1:nLayers))*iden_water
   
 ! end associate statements
 end associate globalVars

 end subroutine opSplittin


 ! **********************************************************************************************************
 ! private subroutine stateSubset: get a mask for the desired state variables
 ! **********************************************************************************************************
 subroutine stateFilter(ixSolution,iStateTypeSplit,iDomainSplit,indx_data,stateMask,nSubset,err,message)
 implicit none
 ! input
 integer(i4b),intent(in)         :: ixSolution                    ! index of solution method (1,2,3,...)
 integer(i4b),intent(in)         :: iStateTypeSplit               ! index of the state type split
 integer(i4b),intent(in)         :: iDomainSplit                  ! index of the domain split
 type(var_ilength),intent(inout) :: indx_data                     ! indices for a local HRU
 ! output
 logical(lgt),intent(out)        :: stateMask(:)                  ! mask defining desired state variables
 integer(i4b),intent(out)        :: nSubset                       ! number of selected state variables for a given split
 integer(i4b),intent(out)        :: err                           ! error code
 character(*),intent(out)        :: message                       ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! data structures
 associate(&
 ! indices of model state variables
 ixStateType => indx_data%var(iLookINDEX%ixStateType)%dat,& ! intent(in): [i4b(:)] indices defining the type of the state (ixNrgState...)
 ixNrgCanair => indx_data%var(iLookINDEX%ixNrgCanair)%dat,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in canopy air space domain
 ixNrgCanopy => indx_data%var(iLookINDEX%ixNrgCanopy)%dat,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the canopy domain
 ixHydCanopy => indx_data%var(iLookINDEX%ixHydCanopy)%dat,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
 ixNrgLayer  => indx_data%var(iLookINDEX%ixNrgLayer)%dat ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer  => indx_data%var(iLookINDEX%ixHydLayer)%dat ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ! number of layers
 nSnow       => indx_data%var(iLookINDEX%nSnow)%dat(1)   ,& ! intent(in): [i4b]    number of snow layers
 nSoil       => indx_data%var(iLookINDEX%nSoil)%dat(1)   ,& ! intent(in): [i4b]    number of soil layers
 nLayers     => indx_data%var(iLookINDEX%nLayers)%dat(1)  & ! intent(in): [i4b]    total number of layers
 ) ! data structures
 ! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='stateFilter/'

 ! identify splitting option
 select case(ixSolution)
 
  ! -----
  ! - fully coupled...
  ! ------------------

  ! use all state variables
  case(fullyCoupled); stateMask(:) = .true.
 
  ! -----
  ! - splitting by state type...
  ! ----------------------------
 
  ! split into energy and mass
  case(splitStateType)
   select case(iStateTypeSplit)
    case(nrgSplit);  stateMask = (ixStateType==iname_nrgCanair .or. ixStateType==iname_nrgCanopy .or. ixStateType==iname_nrgLayer)
    case(massSplit); stateMask = (ixStateType==iname_liqCanopy .or. ixStateType==iname_liqLayer  .or. ixStateType==iname_lmpLayer)
    case default; err=20; message=trim(message)//'unable to identify split based on state type'; return
   end select
 
  ! -----
  ! - splitting by domain...
  ! ------------------------
 
  ! split into vegetation, snow, and soil
  case(splitDomainType,explicitEuler)
 
   ! define state mask
   stateMask=.false. ! (initialize state mask)
   select case(iStateTypeSplit)
 
    ! define mask for energy
    case(nrgSplit)
     select case(iDomainSplit)
      case(vegSplit)
       if(ixNrgCanair(1)/=integerMissing) stateMask(ixNrgCanair) = .true.  ! energy of the canopy air space
       if(ixNrgCanopy(1)/=integerMissing) stateMask(ixNrgCanopy) = .true.  ! energy of the vegetation canopy
       stateMask(ixNrgLayer(1)) = .true.  ! energy of the upper-most layer in the snow+soil domain
      case(snowSplit); if(nSnow>1) stateMask(ixNrgLayer(2:nSnow)) = .true.    ! NOTE: (2:) top layer in the snow+soil domain included in vegSplit
      case(soilSplit); stateMask(ixNrgLayer(max(2,nSnow+1):nLayers)) = .true. ! NOTE: max(2,nSnow+1) gives second layer unless more than 2 snow layers
      case default; err=20; message=trim(message)//'unable to identify model sub-domain'; return
     end select       
 
    ! define mask for water
    case(massSplit) 
     select case(iDomainSplit)
      case(vegSplit);  if(ixHydCanopy(1)/=integerMissing) stateMask(ixHydCanopy) = .true.  ! hydrology of the vegetation canopy
      case(snowSplit); stateMask(ixHydLayer(1:nSnow)) = .true.  ! snow hydrology
      case(soilSplit); stateMask(ixHydLayer(nSnow+1:nLayers)) = .true.  ! soil hydrology
      case default; err=20; message=trim(message)//'unable to identify model sub-domain'; return
     end select
 
    ! check
    case default; err=20; message=trim(message)//'unable to identify the state type'; return
   end select  ! (split based on state type)
 
  case default; err=20; message=trim(message)//'unable to identify solution method'; return
 end select  ! (selecting solution method)
 
 ! get the number of selected state variables
 nSubset = count(stateMask)
 
 ! end associations
 end associate

 end subroutine stateFilter

end module opSplittin_module
