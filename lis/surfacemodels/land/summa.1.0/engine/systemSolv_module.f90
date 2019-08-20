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

module systemSolv_module

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
                    LH_fus,       & ! latent heat of fusion                (J K-1)
                    Tfreeze,      & ! temperature at freezing              (K)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

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
public::systemSolv

! control parameters
real(dp),parameter  :: valueMissing=-9999._dp     ! missing value
real(dp),parameter  :: verySmall=1.e-12_dp        ! a very small number (used to check consistency)
real(dp),parameter  :: veryBig=1.e+20_dp          ! a very big number
real(dp),parameter  :: dx = 1.e-8_dp              ! finite difference increment

contains


 ! **********************************************************************************************************
 ! public subroutine systemSolv: run the coupled energy-mass model for one timestep
 ! **********************************************************************************************************
 subroutine systemSolv(&
                       ! input: model control
                       dt,                & ! intent(in):    time step (s)
                       nState,            & ! intent(in):    total number of state variables
                       firstSubStep,      & ! intent(in):    flag to denote first sub-step
                       firstFluxCall,     & ! intent(inout): flag to indicate if we are processing the first flux call
                       explicitEuler,     & ! intent(in):    flag to denote computing the explicit Euler solution
                       computeVegFlux,    & ! intent(in):    flag to denote if computing energy flux over vegetation
                       ! input/output: data structures
                       type_data,         & ! intent(in):    type of vegetation and soil
                       attr_data,         & ! intent(in):    spatial attributes
                       forc_data,         & ! intent(in):    model forcing data
                       mpar_data,         & ! intent(in):    model parameters
                       indx_data,         & ! intent(inout): index data
                       prog_data,         & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,         & ! intent(inout): model diagnostic variables for a local HRU
                       flux_temp,         & ! intent(inout): model fluxes for a local HRU
                       bvar_data,         & ! intent(in):    model variables for the local basin
                       model_decisions,   & ! intent(in):    model decisions
                       stateVecInit,      & ! intent(in):    initial state vector
                       ! output
                       deriv_data,        & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ixSaturation,      & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       untappedMelt,      & ! intent(out):   un-tapped melt energy (J m-3 s-1)
                       stateVecTrial,     & ! intent(out):   updated state vector
                       explicitError,     & ! intent(out):   error in the explicit solution
                       reduceCoupledStep, & ! intent(out):   flag to reduce the length of the coupled step
                       tooMuchMelt,       & ! intent(out):   flag to denote that there was too much melt 
                       niter,             & ! intent(out):   number of iterations taken
                       err,message)         ! intent(out):   error code and error message
 ! ---------------------------------------------------------------------------------------
 ! structure allocations
 USE globalData,only:flux_meta                        ! metadata on the model fluxes
 USE allocspace_module,only:allocLocal                ! allocate local data structures
 ! simulation of fluxes and residuals given a trial state vector
 USE eval8summa_module,only:eval8summa                ! simulation of fluxes and residuals given a trial state vector
 USE summaSolve_module,only:summaSolve                ! calculate the iteration increment, evaluate the new state, and refine if necessary
 USE getVectorz_module,only:getScaling                ! get the scaling vectors
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                            ! time step (seconds)
 integer(i4b),intent(in)         :: nState                        ! total number of state variables
 logical(lgt),intent(in)         :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall                 ! flag to define the first flux call
 logical(lgt),intent(in)         :: explicitEuler                 ! flag to denote computing the explicit Euler solution
 logical(lgt),intent(in)         :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input/output: data structures
 type(var_i),intent(in)          :: type_data                     ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                     ! spatial attributes
 type(var_d),intent(in)          :: forc_data                     ! model forcing data
 type(var_dlength),intent(in)    :: mpar_data                     ! model parameters
 type(var_ilength),intent(inout) :: indx_data                     ! indices for a local HRU
 type(var_dlength),intent(inout) :: prog_data                     ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                     ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_temp                     ! model fluxes for a local HRU
 type(var_dlength),intent(in)    :: bvar_data                     ! model variables for the local basin
 type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
 real(dp),intent(in)             :: stateVecInit(:)               ! initial state vector (mixed units)
 ! output: model control
 type(var_dlength),intent(inout) :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables 
 integer(i4b),intent(inout)      :: ixSaturation                  ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(dp),intent(out)            :: untappedMelt(:)               ! un-tapped melt energy (J m-3 s-1)
 real(dp),intent(out)            :: stateVecTrial(:)              ! trial state vector (mixed units)
 real(dp),intent(out)            :: explicitError                 ! error in the explicit solution
 logical(lgt),intent(out)        :: reduceCoupledStep             ! flag to reduce the length of the coupled step 
 logical(lgt),intent(out)        :: tooMuchMelt                   ! flag to denote that there was too much melt 
 integer(i4b),intent(out)        :: niter                         ! number of iterations taken
 integer(i4b),intent(out)        :: err                           ! error code
 character(*),intent(out)        :: message                       ! error message
 ! *********************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************
 ! ---------------------------------------------------------------------------------------
 ! * general local variables
 ! ---------------------------------------------------------------------------------------
 character(LEN=256)              :: cmessage                      ! error message of downwind routine
 integer(i4b)                    :: iter                          ! iteration index
 integer(i4b)                    :: iVar                          ! index of variable
 integer(i4b)                    :: iLayer                        ! index of layer in the snow+soil domain
 integer(i4b)                    :: iState                        ! index of model state
 integer(i4b)                    :: nLeadDim                      ! length of the leading dimension of the Jacobian matrix (nBands or nState)
 integer(i4b)                    :: local_ixGroundwater           ! local index for groundwater representation
 real(dp),parameter              :: tempAccelerate=0.00_dp        ! factor to force initial canopy temperatures to be close to air temperature
 real(dp),parameter              :: xMinCanopyWater=0.0001_dp     ! minimum value to initialize canopy water (kg m-2)
 real(dp),parameter              :: tinyStep=0.000001_dp          ! stupidly small time step (s) 
 ! ------------------------------------------------------------------------------------------------------
 ! * model solver
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: forceFullMatrix=.false.       ! flag to force the use of the full Jacobian matrix
 integer(i4b)                    :: maxiter                       ! maximum number of iterations
 integer(i4b)                    :: ixMatrix                      ! form of matrix (band diagonal or full matrix)
 type(var_dlength)               :: flux_init                     ! model fluxes at the start of the time step 
 real(dp),allocatable            :: dBaseflow_dMatric(:,:)        ! derivative in baseflow w.r.t. matric head (s-1)  ! NOTE: allocatable, since not always needed
 real(dp)                        :: stateVecNew(nState)           ! new state vector (mixed units)
 real(dp)                        :: rhsFlux0(nState)              ! right-hand-side fluxes (start of step)
 real(dp)                        :: rhsFlux1(nState)              ! right-hand-side fluxes (end of step)
 real(dp)                        :: cf0(nState),cf1(nState)       ! conversion factor: factor to convert fluxes to states (different flux evaluations)
 real(dp)                        :: fluxVec0(nState)              ! flux vector (mixed units)
 real(dp)                        :: fScale(nState)                ! characteristic scale of the function evaluations (mixed units)
 real(dp)                        :: xScale(nState)                ! characteristic scale of the state vector (mixed units)
 real(dp)                        :: dMat(nState)                  ! diagonal matrix (excludes flux derivatives)
 real(qp)                        :: sMul(nState)    ! NOTE: qp    ! multiplier for state vector for the residual calculations
 real(qp)                        :: rVec(nState)    ! NOTE: qp    ! residual vector
 real(dp)                        :: rAdd(nState)                  ! additional terms in the residual vector
 real(dp)                        :: fOld,fNew                     ! function values (-); NOTE: dimensionless because scaled
 logical(lgt)                    :: stateConstrained              ! flag to denote if the state was constrained in the explicit update
 logical(lgt)                    :: feasible                      ! flag to define the feasibility of the solution
 logical(lgt)                    :: converged                     ! convergence flag
 real(dp)                        :: resSinkNew(nState)            ! additional terms in the residual vector
 real(dp)                        :: fluxVecNew(nState)            ! new flux vector
 real(qp)                        :: resVecNew(nState)  ! NOTE: qp ! new residual vector
 real(dp)                        :: solutionError(nState)         ! vector of errors in the model solution
 real(dp),dimension(1)           :: errorTemp                     ! maximum error in explicit solution
 real(dp)                        :: stateVecUpdate(nState)        ! state vector update
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 globalVars: associate(&
 ! model decisions
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,& ! intent(in):    [i4b]    groundwater parameterization
 ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,& ! intent(in):    [i4b]    spatial representation of groundwater (local-column or single-basin)
 ! accelerate solutuion for temperature
 airtemp                 => forc_data%var(iLookFORCE%airtemp)                      ,& ! intent(in):    [dp]     temperature of the upper boundary of the snow and soil domains (K)
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy energy state variable
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)
 ! vector of energy and hydrology indices for the snow and soil domains
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):    [i4b]    number of hydrology state variables in the snow+soil domain
 ! type of state and domain for a given variable
 ixStateType_subset      => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in):    [i4b(:)] [state subset] type of desired model state variables
 ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat      ,& ! intent(in):    [i4b(:)] [state subset] domain for desired model state variables
 ! layer geometry
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):    [i4b]    number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):    [i4b]    number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in):    [i4b]    total number of layers
 )
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="systemSolv/"

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! -----
 ! * initialize...
 ! ---------------

 ! check
 if(dt < tinyStep)then
  message=trim(message)//'dt is tiny'
  err=20; return
 endif

 ! initialize the flags
 tooMuchMelt        = .false.   ! too much melt
 reduceCoupledStep  = .false.   ! need to reduce the length of the coupled step 

 ! define maximum number of iterations
 maxiter = nint(mpar_data%var(iLookPARAM%maxiter)%dat(1))

 ! modify the groundwater representation for this single-column implementation
 select case(ixSpatialGroundwater)
  case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation)

 ! allocate space for the model fluxes at the start of the time step
 call allocLocal(flux_meta(:),flux_init,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the baseflow derivatives
 ! NOTE: needs allocation because only used when baseflow sinks are active
 if(ixGroundwater==qbaseTopmodel)then
  allocate(dBaseflow_dMatric(nSoil,nSoil),stat=err)  ! baseflow depends on total storage in the soil column, hence on matric head in every soil layer
 else
  allocate(dBaseflow_dMatric(0,0),stat=err)          ! allocate zero-length dimnensions to avoid passing around an unallocated matrix
 end if 
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the baseflow derivatives'; return; end if

 ! identify the matrix solution method
 ! (the type of matrix used to solve the linear system A.X=B)
 if(local_ixGroundwater==qbaseTopmodel .or. forceFullMatrix)then
  nLeadDim=nState         ! length of the leading dimension
  ixMatrix=ixFullMatrix   ! named variable to denote the full Jacobian matrix
 else
  nLeadDim=nBands         ! length of the leading dimension
  ixMatrix=ixBandMatrix   ! named variable to denote the band-diagonal matrix
 endif
   
 ! initialize the model fluxes (some model fluxes are not computed in the iterations)
 do iVar=1,size(flux_temp%var)
  flux_init%var(iVar)%dat(:) = flux_temp%var(iVar)%dat(:)
 end do
   
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 ! *** NUMERICAL SOLUTION FOR A GIVEN SUBSTEP AND SPLIT *********************************************************************
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 ! **************************************************************************************************************************
 
 ! -----
 ! * get scaling vectors...
 ! ------------------------
 
 ! initialize state vectors
 call getScaling(&
                 ! input
                 diag_data,                        & ! intent(in):    model diagnostic variables for a local HRU
                 indx_data,                        & ! intent(in):    indices defining model states and layers
                 ! output
                 fScale,                           & ! intent(out):   function scaling vector (mixed units)
                 xScale,                           & ! intent(out):   variable scaling vector (mixed units)
                 sMul,                             & ! intent(out):   multiplier for state vector (used in the residual calculations)
                 dMat,                             & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes) 
                 err,cmessage)                       ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
 
 ! -----
 ! * compute the initial function evaluation...
 ! --------------------------------------------
 
 ! initialize the trial state vectors
 stateVecTrial = stateVecInit
 
 ! need to intialize canopy water at a positive value
 if(ixVegHyd/=integerMissing)then
  if(stateVecTrial(ixVegHyd) < xMinCanopyWater) stateVecTrial(ixVegHyd) = stateVecTrial(ixVegHyd) + xMinCanopyWater
 endif
 
 ! try to accelerate solution for energy
 if(ixCasNrg/=integerMissing) stateVecTrial(ixCasNrg) = stateVecInit(ixCasNrg) + (airtemp - stateVecInit(ixCasNrg))*tempAccelerate
 if(ixVegNrg/=integerMissing) stateVecTrial(ixVegNrg) = stateVecInit(ixVegNrg) + (airtemp - stateVecInit(ixVegNrg))*tempAccelerate

 ! compute the flux and the residual vector for a given state vector
 ! NOTE 1: The derivatives computed in eval8summa are used to calculate the Jacobian matrix for the first iteration
 ! NOTE 2: The Jacobian matrix together with the residual vector is used to calculate the first iteration increment
 call eval8summa(&
                 ! input: model control
                 dt,                      & ! intent(in):    length of the time step (seconds)
                 nSnow,                   & ! intent(in):    number of snow layers
                 nSoil,                   & ! intent(in):    number of soil layers
                 nLayers,                 & ! intent(in):    number of layers
                 nState,                  & ! intent(in):    number of state variables in the current subset
                 firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                 firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                 .true.,                  & ! intent(in):    flag to indicate if we are processing the first iteration in a splitting operation
                 computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                 ! input: state vectors
                 stateVecTrial,           & ! intent(in):    model state vector
                 fScale,                  & ! intent(in):    function scaling vector
                 sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
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
                 flux_init,               & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                 deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! input-output: baseflow
                 ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                 dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                 ! output
                 feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                 fluxVec0,                & ! intent(out):   flux vector
                 rAdd,                    & ! intent(out):   additional (sink) terms on the RHS of the state equation
                 rVec,                    & ! intent(out):   residual vector
                 fOld,                    & ! intent(out):   function evaluation
                 err,cmessage)              ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! check feasibility (state vector SHOULD be feasible at this point)
 if(.not.feasible)then
  reduceCoupledStep=.true.
  return
 endif
 
 ! copy over the initial flux structure since some model fluxes are not computed in the iterations
 do concurrent ( iVar=1:size(flux_meta) )
  flux_temp%var(iVar)%dat(:) = flux_init%var(iVar)%dat(:)
 end do
 
 ! ** if explicit Euler, then estimate state vector at the end of the time step
 if(explicitEuler)then

  ! --> compute the RHS fluxes and conversion factor
  call rhsFluxes(indx_data,deriv_data,sMul,fluxVec0,rAdd/dt,      & ! intent(in)  : state indices and derivatives, and the state vector multiplier
                 cf0,rhsFlux0,err,cmessage)                         ! intent(out) : conversion factor, and error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! --> update states using the explicit Euler method
  call explicitUpdate(indx_data,mpar_data,prog_data,stateVecInit, & ! intent(in)  : indices, parameters, prognostic variables, and initial state vector
                      dt*fluxVec0/cf0,                            & ! intent(in)  : state vector update      
                      stateVecTrial,stateConstrained,             & ! intent(out) : trial state vector and flag to denote that it was constrained
                      err,cmessage)                                 ! intent(out) : error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 endif  ! if explicit Euler
 
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 
 ! **************************
 ! *** MAIN ITERATION LOOP...
 ! **************************
 
 ! iterate
 ! NOTE: this do loop is skipped in the explicitEuler solution (localMaxIter=0)
 do iter=1,maxIter
 
  ! print iteration count
  !print*, '*** iter, maxiter, dt = ', iter, maxiter, dt

  ! keep track of the number of iterations
  niter = iter+1  ! +1 because xFluxResid was moved outside the iteration loop (for backwards compatibility)

  ! compute the next trial state vector
  !  1) Computes the Jacobian matrix based on derivatives from the last flux evaluation
  !  2) Computes the iteration increment based on Jacobian and residuals from the last flux evaluation
  !  3) Computes new fluxes and derivatives, new residuals, and (if necessary) refines the state vector
  ! NOTE: only returns the flux vector and function evaluation when the solution method is explicitEuler
  call summaSolve(&
                  ! input: model control
                  dt,                            & ! intent(in):    length of the time step (seconds)
                  explicitEuler,                 & ! intent(in):    logical flag to only return the flux and function evaluation
                  iter,                          & ! intent(in):    iteration index
                  nSnow,                         & ! intent(in):    number of snow layers
                  nSoil,                         & ! intent(in):    number of soil layers
                  nLayers,                       & ! intent(in):    total number of layers
                  nLeadDim,                      & ! intent(in):    length of the leading dimension of the Jacobian matrix (either nBands or nState) 
                  nState,                        & ! intent(in):    total number of state variables
                  ixMatrix,                      & ! intent(in):    type of matrix (full or band diagonal)
                  firstSubStep,                  & ! intent(in):    flag to indicate if we are processing the first sub-step
                  firstFluxCall,                 & ! intent(inout): flag to indicate if we are processing the first flux call
                  computeVegFlux,                & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  ! input: state vectors
                  stateVecTrial,                 & ! intent(in):    trial state vector
                  fScale,                        & ! intent(in):    function scaling vector
                  xScale,                        & ! intent(in):    "variable" scaling vector, i.e., for state variables
                  rVec,                          & ! intent(in):    residual vector
                  sMul,                          & ! intent(in):    state vector multiplier (used in the residual calculations)
                  dMat,                          & ! intent(inout): diagonal matrix (excludes flux derivatives)
                  fOld,                          & ! intent(in):    old function evaluation
                  ! input: data structures       
                  model_decisions,               & ! intent(in):    model decisions
                  type_data,                     & ! intent(in):    type of vegetation and soil
                  attr_data,                     & ! intent(in):    spatial attributes
                  mpar_data,                     & ! intent(in):    model parameters
                  forc_data,                     & ! intent(in):    model forcing data
                  bvar_data,                     & ! intent(in):    average model variables for the entire basin
                  prog_data,                     & ! intent(in):    model prognostic variables for a local HRU
                  indx_data,                     & ! intent(in):    index data
                  ! input-output: data structures
                  diag_data,                     & ! intent(inout): model diagnostic variables for a local HRU
                  flux_temp,                     & ! intent(inout): model fluxes for a local HRU (temporary structure)
                  deriv_data,                    & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                  ! input-output: baseflow       
                  ixSaturation,                  & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                  dBaseflow_dMatric,             & ! intent(inout): derivative in baseflow w.r.t. matric head (s-1)
                  ! output
                  stateVecNew,                   & ! intent(out):   new state vector
                  fluxVecNew,                    & ! intent(out):   new flux vector
                  resSinkNew,                    & ! intent(out):   additional (sink) terms on the RHS of the state equa
                  resVecNew,                     & ! intent(out):   new residual vector
                  fNew,                          & ! intent(out):   new function evaluation
                  feasible,                      & ! intent(out):   flag to denote that the state vector is feasible
                  converged,                     & ! intent(out):   convergence flag
                  err,cmessage)                    ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  !print*, err,trim(cmessage) 

  ! update function evaluation, residual vector, and states
  ! NOTE 1: The derivatives computed in summaSolve are used to calculate the Jacobian matrix at the next iteration
  ! NOTE 2: The Jacobian matrix together with the residual vector is used to calculate the new iteration increment
  if(.not.explicitEuler)then
   ! save functions and residuals
   fOld          = fNew
   rVec          = resVecNew 
   stateVecTrial = stateVecNew
   ! check feasibility
   if(.not.feasible)then
    message=trim(message)//'expect feasible solution in implicit Euler'
    err=20; return
   endif  ! check feasibility
  endif  ! check explicit Euler
 
  ! print progress
  !write(*,'(a,10(f16.14,1x))') 'rVec                  = ', rVec           ( min(nState,iJac1) : min(nState,iJac2) )
  !write(*,'(a,10(f16.10,1x))') 'fluxVecNew            = ', fluxVecNew     ( min(nState,iJac1) : min(nState,iJac2) )*dt
  !write(*,'(a,10(f16.10,1x))') 'stateVecTrial         = ', stateVecTrial  ( min(nState,iJac1) : min(nState,iJac2) )
  !print*, 'PAUSE: check states and fluxes'; read(*,*) 
 
  ! exit iteration loop if converged
  if(converged .or. explicitEuler) exit
 
  ! check convergence
  if(iter==maxiter)then
   message=trim(message)//'failed to converge'
   err=-20; return
  endif
  !print*, 'PAUSE: iterating'; read(*,*)
 
 end do  ! iterating
 !print*, 'PAUSE: after iterations'; read(*,*)
  
 ! -----
 ! * update states...
 ! ------------------
 
 ! special case of explicit Euler
 if(explicitEuler)then

  ! --> check feasibility
  if(.not.feasible)then
   message=trim(message)//'state is not feasible in explicit Euler (reduce time step)'
   err=-20; return
  endif

  ! --> compute the RHS fluxes and conversion factor
  call rhsFluxes(indx_data,deriv_data,sMul,fluxVecNew,resSinkNew/dt, &  ! intent(in)  : state indices and derivatives, and the state vector multiplier
                 cf1,rhsFlux1,err,cmessage)                             ! intent(out) : conversion factor, and error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! --> compute state vector update
  stateVecUpdate = 0.5_dp*( dt*fluxVec0/cf0 + dt*fluxVecNew/cf1 )

  ! compute melt energy for the explicit Euler method
  call explicitMelt(dt,indx_data,diag_data,prog_data,     &  ! intent(in)    : time step and data structures 
                    0.5_dp*(fluxVec0 + fluxVecNew),sMul,  &  ! intent(in)    : total flux, and state vector multipler
                    stateVecUpdate,                       &  ! intent(inout) : state vector update (modified if ice cannot support melt)
                    untappedMelt,                         &  ! intent(out)   : untapped melt energy (J m-3 s-1)
                    tooMuchMelt,                          &  ! intent(out)   : flag to denote that ice is insufficient to support available melt
                    err,cmessage)                            ! intent(out)   : error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! if too much melt then return
  ! NOTE: need to go all the way back to coupled_em and merge snow layers, as all splitting operations need to occur with the same layer geometry
  if(tooMuchMelt)then
   reduceCoupledStep=.true.
   return
  endif

  ! --> update states using the explicit Euler method
  call explicitUpdate(indx_data,mpar_data,prog_data,            &  ! intent(in)  : indices, parameters, prognostic variables
                      stateVecInit,stateVecUpdate,stateVecTrial,&  ! intent(in)  : initial state vector, state update, and new state vector
                      stateConstrained,err,cmessage)               ! intent(out) : flag for state contraint, and error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! --> estimate the solution error
  ! NOTE: done before the constraints check to return the error
  solutionError(:) = abs(fluxVec0*dt - cf1*stateVecUpdate)
  errorTemp        = maxval(solutionError)
  explicitError    = max(errorTemp(1), verySmall)

  ! print progress in the explicit Euler solution
  if(globalPrintFlag)then
   write(*,'(a,1x,10(f20.12,1x))') 'cf0           = ', cf0            ( min(nState,iJac1) : min(nState,iJac2) )
   write(*,'(a,1x,10(f20.12,1x))') 'cf1           = ', cf1            ( min(nState,iJac1) : min(nState,iJac2) )
   write(*,'(a,1x,10(f20.12,1x))') 'fluxVec0      = ', fluxVec0       ( min(nState,iJac1) : min(nState,iJac2) )
   write(*,'(a,1x,10(f20.12,1x))') 'fluxVecNew    = ', fluxVecNew     ( min(nState,iJac1) : min(nState,iJac2) )
   write(*,'(a,1x,10(f20.12,1x))') 'rAdd          = ', rAdd           ( min(nState,iJac1) : min(nState,iJac2) )
   write(*,'(a,1x,10(f20.12,1x))') 'resSinkNew    = ', resSinkNew     ( min(nState,iJac1) : min(nState,iJac2) )
   write(*,'(a,1x,10(f20.12,1x))') 'stateVecInit  = ', stateVecInit   ( min(nState,iJac1) : min(nState,iJac2) )
   write(*,'(a,1x,10(f20.12,1x))') 'stateVecTrial = ', stateVecTrial  ( min(nState,iJac1) : min(nState,iJac2) )
   write(*,'(a,1x,10(f20.12,1x))') 'stateVecNew   = ', stateVecNew    ( min(nState,iJac1) : min(nState,iJac2) )
   write(*,'(a,1x,10(f20.12,1x))') 'solutionError = ', solutionError  ( min(nState,iJac1) : min(nState,iJac2) )
   print*, 'dt = ', dt
   !print*, 'PAUSE: checking state vector for the explicit Euler solution'; read(*,*)
  endif  ! (if printing)

  

  ! check if the state is constrained
  if(stateConstrained)then
   message=trim(message)//'state is constrained in explicit Heun (reduce time step)'
   err=-20; return
  endif

  ! average start-of-step and end-of-step fluxes for explicit Euler
  do iVar=1,size(flux_temp%var)
   flux_temp%var(iVar)%dat(:) = 0.5_dp*(flux_init%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:) )
  end do
   
 ! standard implicit solution
 else  ! switch between explicit and implicit Euler

  ! set explicit error to missing
  explicitError = realMissing

  ! set untapped melt energy to zero
  untappedMelt(:) = 0._dp

  ! update temperatures (ensure new temperature is consistent with the fluxes)
  if(nSnowSoilNrg>0)then
   do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
    iState = ixSnowSoilNrg(iLayer)
    stateVecTrial(iState) = stateVecInit(iState) + (fluxVecNew(iState)*dt + resSinkNew(iState))/real(sMul(iState), dp)
   end do  ! looping through non-missing energy state variables in the snow+soil domain
  endif

  ! update volumetric water content in the snow (ensure change in state is consistent with the fluxes)
  ! NOTE: for soil water balance is constrained within the iteration loop
  if(nSnowSoilHyd>0)then
   do concurrent (iLayer=1:nSnow,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow domain)
    iState = ixSnowSoilHyd(iLayer)
    stateVecTrial(iState) = stateVecInit(iState) + (fluxVecNew(iState)*dt + resSinkNew(iState))
   end do  ! looping through non-missing energy state variables in the snow+soil domain
  endif

 endif  ! switch between explicit and implicit Euler
      
 ! end associate statements
 end associate globalVars

 end subroutine systemSolv

 ! **********************************************************************************************************
 ! private subroutine: compute the right-hand-side fluxes and conversion factors 
 ! **********************************************************************************************************
 subroutine rhsFluxes(&
                      ! input
                      indx_data,                    & ! intent(in)  : state indices
                      deriv_data,                   & ! intent(in)  : state derivatives
                      stateVecMult,                 & ! intent(in)  : state vector multiplier
                      fluxVec,                      & ! intent(in)  : flux vector
                      sink,                         & ! intent(in)  : sink
                      ! output
                      conversionFactor,             & ! intent(out) : flux2state conversion factor
                      rhsFlux,                      & ! intent(out) : right-hand-side fluxes
                      err,message)                    ! intent(out) : error control
 USE var_lookup,only:iLookINDEX                       ! named variables for structure elements
 USE var_lookup,only:iLookDERIV                       ! named variables for structure elements
 implicit none
 ! input
 type(var_ilength),intent(in)  :: indx_data           ! state indices
 type(var_dlength),intent(in)  :: deriv_data          ! derivatives in model fluxes w.r.t. relevant state variables
 real(qp)         ,intent(in)  :: stateVecMult(:)     ! state vector multiplier 
 real(dp)         ,intent(in)  :: fluxVec(:)          ! flux vector
 real(dp)         ,intent(in)  :: sink(:)             ! sink
 ! output
 real(dp)         ,intent(out) :: conversionFactor(:) ! change in state w.r.t. time (dS/dt)
 real(dp)         ,intent(out) :: rhsFlux(:)          ! right-hand-side flux
 integer(i4b)     ,intent(out) :: err                 ! error code
 character(*)     ,intent(out) :: message             ! error message
 ! local variables
 integer(i4b)                  :: iState              ! state index
 integer(i4b)                  :: ixFullVector        ! index in the full state vector
 integer(i4b)                  :: ixControlIndex      ! index of the control volume for different domains (veg, snow, soil) 
 real(dp)                      :: meltDeriv           ! melt derivative (J m-3 K-1)
 ! make association with model indices and model derivatives
 associate(&
  ! derivatives
  dVolTot_dPsi0       => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat,       & ! intent(in): [dp(:)]  derivative in total water content w.r.t. total water matric potential
  dTheta_dTkCanopy    => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1), & ! intent(in): [dp]     derivative of volumetric liquid water content w.r.t. temperature (K-1)
  mLayerdTheta_dTk    => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat,    & ! intent(in): [dp(:)]  derivative of volumetric liquid water content w.r.t. temperature (K-1)
  ! indices
  nSnow               => indx_data%var(iLookINDEX%nSnow)%dat(1),             & ! intent(in): [i4b]    number of snow layers
  ixControlVolume     => indx_data%var(iLookINDEX%ixControlVolume)%dat,      & ! intent(in): [i4b(:)] index of the control volume for different domains (veg, snow, soil)
  ixMapSubset2Full    => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat,     & ! intent(in): [i4b(:)] [state subset] list of indices of the full state vector in the state subset
  ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat,   & ! intent(in): [i4b(:)] [state subset] type of desired model state variables
  ixDomainType_subset => indx_data%var(iLookINDEX%ixDomainType_subset)%dat   & ! intent(in): [i4b(:)] [state subset] type of desired model state variables
 ) ! associations
 ! ---------------------------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='rhsFluxes/'

 ! loop through model states
 do iState=1,size(stateVecMult)

  ! get index of the control volume within the domain
  ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
  ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

  ! if ice included in the energy equation, then get the melt derivative
  if(ixStateType_subset(iState)==iname_nrgCanopy .or. ixStateType_subset(iState)==iname_nrgLayer)then
   select case( ixDomainType_subset(iState) )
    case(iname_veg);  meltDeriv = LH_fus*iden_water*dTheta_dTkCanopy 
    case(iname_snow); meltDeriv = LH_fus*iden_water*mLayerdTheta_dTk(ixControlIndex)
    case(iname_soil); meltDeriv = LH_fus*iden_water*mLayerdTheta_dTk(ixControlIndex+nSnow)
    case default; err=20; message=trim(message)//'expect domain type to be iname_veg, iname_snow or iname_soil'; return
   end select
  else
   meltDeriv = 0._dp
  endif

  ! get the flux-2-state conversion factor
  select case( ixStateType_subset(iState) )
   case(iname_matLayer,iname_lmpLayer);                                 conversionFactor(iState) = dVolTot_dPsi0(ixControlIndex)
   case(iname_nrgCanair,iname_nrgCanopy,iname_nrgLayer);                conversionFactor(iState) = real(stateVecMult(iState), dp) + meltDeriv
   case(iname_watCanopy,iname_liqCanopy,iname_watLayer,iname_liqLayer); conversionFactor(iState) = 1._dp
   case default; err=20; message=trim(message)//'unable to identify the state type'; return
  end select

  ! get the right-hand-side flux
  select case( ixStateType_subset(iState) )
   case(iname_matLayer,iname_lmpLayer,iname_watLayer,iname_liqLayer);   rhsFlux(iState) = fluxVec(iState) + sink(iState) ! include transpiration and lateral flow sinks as part of the flux
   case(iname_nrgCanair,iname_nrgCanopy,iname_nrgLayer);                rhsFlux(iState) = fluxVec(iState) 
   case(iname_watCanopy,iname_liqCanopy);                               rhsFlux(iState) = fluxVec(iState) 
   case default; err=20; message=trim(message)//'unable to identify the state type'; return
  end select

 end do  ! looping through state variables

 ! end association with data structures
 end associate

 end subroutine rhsFluxes

 ! **********************************************************************************************************
 ! private subroutine explicitMelt: calaculate the melt in the explicit solution
 ! **********************************************************************************************************
 subroutine explicitMelt(&
                         dt,                        & ! intent(in)    : time step (s)
                         indx_data,                 & ! intent(in)    : state indices
                         diag_data,                 & ! intent(in)    : model diagnostic variables
                         prog_data,                 & ! intent(in)    : model prognostic variables
                         totalFlux,                 & ! intent(in)    : total flux
                         stateVecMult,              & ! intent(in)    : state vector multiplier
                         stateVecUpdate,            & ! intent(inout) : state vector update
                         untappedMelt,              & ! intent(out)   : untapped melt energy (J m-3 s-1)
                         tooMuchMelt,               & ! intent(out)   : flag to denote that ice is insufficient to support available melt
                         err,message)                 ! intent(out)   : error control
 USE var_lookup,only:iLookDIAG                        ! named variables for structure elements
 USE var_lookup,only:iLookPROG                        ! named variables for structure elements
 USE var_lookup,only:iLookINDEX                       ! named variables for structure elements
 implicit none
 ! input
 real(dp)         ,intent(in)    :: dt                ! time step (s)
 type(var_ilength),intent(in)    :: indx_data         ! state indices
 type(var_dlength),intent(in)    :: diag_data         ! model diagnostic variables
 type(var_dlength),intent(in)    :: prog_data         ! model prognostic variables
 real(dp)         ,intent(in)    :: totalFlux(:)      ! total flux
 real(qp)         ,intent(in)    :: stateVecMult(:)   ! state vector multiplier 
 ! output
 real(dp)         ,intent(inout) :: stateVecUpdate(:) ! state vector update
 real(dp)         ,intent(out)   :: untappedMelt(:)   ! untapped melt energy (J m-3 s-1)
 logical(lgt)     ,intent(out)   :: tooMuchMelt       ! flag to denote that ice is insufficient to support available melt
 integer(i4b)     ,intent(out)   :: err               ! error code
 character(*)     ,intent(out)   :: message           ! error message
 ! local variables
 integer(i4b)                    :: iState            ! state index
 integer(i4b)                    :: ixFullVector      ! index in the full state vector
 integer(i4b)                    :: ixControlIndex    ! index of the control volume for different domains (veg, snow, soil)
 real(dp)                        :: tempNrg           ! energy associated with the temperature increase (J m-3 s-1)
 real(dp)                        :: nrg2meltIce       ! energy required to melt all of the ice (J m-3)
 real(dp)                        :: nrg2freezeWater   ! energy required to freeze all of the liquid water (J m-3)
 real(dp)                        :: untappedNrg       ! untapped energy (J m-3)
 real(dp)                        :: xIce              ! ice at the start of the step (kg m-2 [canopy] or dimensionless [snow+soil])
 real(dp)                        :: xLiq              ! liquid water at the start of the step (kg m-2 [canopy] or dimensionless [snow+soil])
 ! --------------------------------------------------------------------------------------------------------------
 ! make association with model indices defined in indexSplit
 associate(&
  ! diagnostic and prognostic variables
  canopyDepth         => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1),  & ! intent(in): [dp]     canopy depth (m)
  scalarCanopyIce     => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),    & ! intent(in): [dp]     ice stored on the vegetation canopy (-)
  scalarCanopyLiq     => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),    & ! intent(in): [dp]     liquid water stored on the vegetation canopy (-)
  mLayerVolFracIce    => prog_data%var(iLookPROG%mLayerVolFracIce)%dat,      & ! intent(in): [dp(:)]  volumetric fraction of ice (-)
  mLayerVolFracLiq    => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat,      & ! intent(in): [dp(:)]  volumetric fraction of liquid water (-)
  ! indices
  nSnow               => indx_data%var(iLookINDEX%nSnow)%dat(1),             & ! intent(in): [i4b]    number of snow layers
  ixControlVolume     => indx_data%var(iLookINDEX%ixControlVolume)%dat,      & ! intent(in): [i4b(:)] index of the control volume for different domains (veg, snow, soil)
  ixMapSubset2Full    => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat,     & ! intent(in): [i4b(:)] [state subset] list of indices of the full state vector in the state subset
  ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat,   & ! intent(in): [i4b(:)] [state subset] type of desired model state variables
  ixDomainType_subset => indx_data%var(iLookINDEX%ixDomainType_subset)%dat   & ! intent(in): [i4b(:)] [state subset] type of desired model state variables
 ) ! associations

 ! initialize error control
 err=0; message='explicitMelt/'

 ! initialize the flag to denote that ice is insufficient to support available melt
 tooMuchMelt=.false.

 ! loop through model states
 do iState=1,size(totalFlux)

  ! --> get index of the control volume within the domain
  ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
  ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

  ! restrict attention to the energy state variables in domains where ice can me be present
  if(ixStateType_subset(iState)==iname_nrgCanopy .or. ixStateType_subset(iState)==iname_nrgLayer)then

   ! --> compute the un-tapped melt energy
   tempNrg              = stateVecUpdate(iState)*real(stateVecMult(iState), dp)/dt  ! energy associated with the temperature increase (J m-3 s-1)
   untappedMelt(iState) = totalFlux(iState) - tempNrg

   ! *****
   ! melting
   if(untappedMelt(iState) > 0._dp)then

    ! --> get the ice at the start of the time step
    select case( ixDomainType_subset(iState) )
     case(iname_veg);  xIce = scalarCanopyIce                             ! kg m-2
     case(iname_snow); xIce = mLayerVolFracIce(ixControlIndex)            ! (-)
     case(iname_soil); xIce = mLayerVolFracIce(ixControlIndex+nSnow)      ! (-)
     case default; err=20; message=trim(message)//'cannot find the domain'; return
    end select
  
    ! --> get the energy required to melt all of the ice (J m-3)
    if(xIce > epsilon(dt))then
     select case( ixDomainType_subset(iState) )
      case(iname_veg);  nrg2meltIce =            LH_fus*xIce/canopyDepth  ! J m-3
      case(iname_snow); nrg2meltIce = iden_ice  *LH_fus*xIce              ! J m-3
      case(iname_soil); nrg2meltIce = iden_water*LH_fus*xIce              ! J m-3
      case default; err=20; message=trim(message)//'cannot find the domain'; return
     end select
    else
     nrg2meltIce = 0._dp
    endif

    ! check if the required melt can be satisfied by the available ice
    if(untappedMelt(iState)*dt > nrg2meltIce)then
  
     ! domain-specfic adjustments
     select case( ixDomainType_subset(iState) )
  
      ! --> vegetation and soil have physical structure, so can recover
      case(iname_veg, iname_soil)
       untappedNrg            = untappedMelt(iState)*dt - nrg2meltIce     ! extra energy not used in melt (J m-3)
       untappedMelt(iState)   = nrg2meltIce/dt                            ! truncate melt to the energy required to melt all ice (J m-3 s-1)
       stateVecUpdate(iState) = stateVecUpdate(iState) + untappedNrg/real(stateVecMult(iState), dp)  ! use the extra energy to update the state vector

      ! --> snow is a problem, as we cannot melt all of the ice in a single time step
      case(iname_snow)
       tooMuchMelt            = .true.
  
      ! --> checks
      case default; err=20; message=trim(message)//'cannot find the domain'; return
     end select

    endif  ! if melt is less than that required to melt all of the ice

   ! *****
   ! freezing
   else

    ! --> get the liquid water at the start of the time step
    select case( ixDomainType_subset(iState) )
     case(iname_veg);  xLiq = scalarCanopyLiq                             ! kg m-2
     case(iname_snow); xLiq = mLayerVolFracLiq(ixControlIndex)            ! (-)
     case(iname_soil); xLiq = mLayerVolFracLiq(ixControlIndex+nSnow)      ! (-)
     case default; err=20; message=trim(message)//'cannot find the domain'; return
    end select

    ! --> get the energy required to freeze all of the liquid water (J m-3)
    if(xLiq > epsilon(dt))then
     select case( ixDomainType_subset(iState) )
      case(iname_veg);  nrg2freezeWater =            LH_fus*xLiq/canopyDepth  ! J m-3
      case(iname_snow); nrg2freezeWater = iden_water*LH_fus*xLiq              ! J m-3
      case(iname_soil); nrg2freezeWater = iden_water*LH_fus*xLiq              ! J m-3
      case default; err=20; message=trim(message)//'cannot find the domain'; return
     end select
    else
     nrg2freezeWater = 0._dp
    endif

    ! check if the required melt can be satisfied by the available ice
    ! NOTE 1: negative untappedMelt
    ! NOTE 2: insufficient liquid water to freeze is never a problem as temperatures just decrease
    if(-untappedMelt(iState)*dt > nrg2freezeWater)then
     untappedNrg            = -untappedMelt(iState)*dt - nrg2freezeWater     ! extra energy not used in melt (J m-3)
     untappedMelt(iState)   = -nrg2freezeWater/dt                            ! truncate melt to the energy required to melt all ice (J m-3 s-1)
     stateVecUpdate(iState) = stateVecUpdate(iState) - untappedNrg/real(stateVecMult(iState), dp)  ! use the extra energy to update the state vector
    endif  ! if freeze is greater than that required to freeze all of the water

   endif    ! if freezing

  ! not a relevant energy state (or not an energy state at all!)
  else
   untappedMelt(iState) = 0._dp
  endif

 end do  ! looping through state variables

 ! end association with data structures
 end associate

 end subroutine explicitMelt


 ! **********************************************************************************************************
 ! private subroutine explicitUpdate: update the states using the explicit Euler method
 ! **********************************************************************************************************
 subroutine explicitUpdate(&
                           indx_data,             & ! intent(in)  : state indices
                           mpar_data,             & ! intent(in)  : model parameters
                           prog_data,             & ! intent(in)  : model prognostic variables
                           stateVecInit,          & ! intent(in)  : initial state vector
                           stateVecUpdate,        & ! intent(in)  : state vector update
                           stateVecNew,           & ! intent(out) : new state vector
                           constrained,           & ! intent(out) : flag to denote if the state was constrained
                           err,message)             ! intent(out) : error control
 USE var_lookup,only:iLookPROG                      ! named variables for structure elements
 USE var_lookup,only:iLookPARAM                     ! named variables for structure elements
 USE var_lookup,only:iLookINDEX                     ! named variables for structure elements
 implicit none
 ! input
 type(var_ilength),intent(in)  :: indx_data         ! state indices
 type(var_dlength),intent(in)  :: mpar_data         ! model parameters
 type(var_dlength),intent(in)  :: prog_data         ! model prognostic variables
 real(dp)         ,intent(in)  :: stateVecInit(:)   ! initial state vector   
 real(dp)         ,intent(in)  :: stateVecUpdate(:) ! state vector update
 ! output
 real(dp)         ,intent(out) :: stateVecNew(:)    ! new state vector
 logical(lgt)     ,intent(out) :: constrained       ! flag to denote if the state was constrained
 integer(i4b)     ,intent(out) :: err               ! error code
 character(*)     ,intent(out) :: message           ! error message
 ! local variables
 integer(i4b)                  :: iState            ! state index
 integer(i4b)                  :: ixFullVector      ! index in the full state vector
 integer(i4b)                  :: ixControlIndex    ! index of the control volume for different domains (veg, snow, soil) 
 real(dp)                      :: valueMin,valueMax ! minimum and maximum state values    
 ! --------------------------------------------------------------------------------------------------------------

 ! make association with model indices defined in indexSplit
 associate(&
  theta_sat           => mpar_data%var(iLookPARAM%theta_sat)%dat,            & ! intent(in): [dp]     soil porosity (-)
  theta_res           => mpar_data%var(iLookPARAM%theta_res)%dat,            & ! intent(in): [dp]     soil residual volumetric water content (-)
  mLayerVolFracIce    => prog_data%var(iLookPROG%mLayerVolFracIce)%dat,      & ! intent(in): [dp(:)]  volumetric fraction of ice (-)
  ixControlVolume     => indx_data%var(iLookINDEX%ixControlVolume)%dat,      & ! intent(in): [i4b(:)] index of the control volume for different domains (veg, snow, soil)
  ixMapSubset2Full    => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat,     & ! intent(in): [i4b(:)] [state subset] list of indices of the full state vector in the state subset
  ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat,   & ! intent(in): [i4b(:)] [state subset] type of desired model state variables
  ixDomainType_subset => indx_data%var(iLookINDEX%ixDomainType_subset)%dat   & ! intent(in): [i4b(:)] [state subset] type of desired model state variables
 ) ! associations

 ! initialize error control
 err=0; message='explicitUpdate/'

 ! initialize the flag to denote if the state is constrained
 constrained=.false.

 ! loop through model states
 do iState=1,size(stateVecInit)

  ! get index of the control volume within the domain
  ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
  ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

  ! update the state vector 
  stateVecNew(iState) = stateVecInit(iState) + stateVecUpdate(iState)

  ! impose non-negativity constraints for the mass of water on the vegetation canopy
  if(ixStateType_subset(iState)==iname_watCanopy .or. ixStateType_subset(iState)==iname_liqCanopy)then
   if(stateVecNew(iState) < 0._dp)then
    stateVecNew(iState)=0._dp
    constrained=.true.
   endif
  endif

  ! impose minimum and maximum storage constraints for volumetric water
  if(ixStateType_subset(iState)==iname_watLayer .or. ixStateType_subset(iState)==iname_liqLayer)then
   select case( ixDomainType_subset(iState) )
    case(iname_snow)
     valueMin = 0._dp
     valueMax = merge(iden_ice/iden_water, 1._dp - mLayerVolFracIce(ixControlIndex), ixStateType_subset(iState)==iname_watLayer)
    case(iname_soil)
     valueMin = theta_res(ixControlIndex)
     valueMax = theta_sat(ixControlIndex)
    case default; err=20; message=trim(message)//'expect domain type to be iname_snow or iname_soil'; return
   end select
   if(stateVecNew(iState) < valueMin)then
    stateVecNew(iState)=valueMin
    constrained=.true.
   endif
   if(stateVecNew(iState) > valueMax)then
    stateVecNew(iState)=valueMax
    constrained=.true.
   endif
  endif

  ! impose below-freezing constraints for snow temperature
  if(ixDomainType_subset(iState)==iname_snow .and. ixStateType_subset(iState)==iname_nrgLayer)then
   if(stateVecNew(iState) > Tfreeze)then
    stateVecNew(iState)=Tfreeze
    constrained=.true.
   endif
  endif

 end do ! looping through states

 ! end association to the information in the data structures
 end associate

 end subroutine explicitUpdate





end module systemSolv_module
