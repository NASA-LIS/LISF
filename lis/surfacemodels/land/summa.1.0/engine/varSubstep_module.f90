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

module varSubstep_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access the global print flag
USE globalData,only:globalPrintFlag

! domain types
USE globalData,only:iname_cas       ! named variables for the canopy air space
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! constants
USE multiconst,only:&
                    Tfreeze,      & ! freezing temperature                 (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! safety: set private unless specified otherwise
implicit none
private
public::varSubstep

! algorithmic parameters
real(dp),parameter     :: verySmall=1.e-6_dp   ! used as an additive constant to check if substantial difference among real numbers

contains


 ! **********************************************************************************************************
 ! public subroutine varSubstep: run the model for a collection of substeps for a given state subset
 ! **********************************************************************************************************
 subroutine varSubstep(&
                       ! input: model control
                       dt,                & ! intent(in)    : time step (s)
                       dt_min,            & ! intent(in)    : minimum time step (seconds)
                       errTol,            & ! intent(in)    : error tolerance for the explicit solution
                       nState,            & ! intent(in)    : total number of state variables
                       doAdjustTemp,      & ! intent(in)    : flag to indicate if we adjust the temperature
                       firstSubStep,      & ! intent(in)    : flag to denote first sub-step
                       firstFluxCall,     & ! intent(inout) : flag to indicate if we are processing the first flux call
                       explicitEuler,     & ! intent(in)    : flag to denote computing the explicit Euler solution
                       computeVegFlux,    & ! intent(in)    : flag to denote if computing energy flux over vegetation
                       fluxMask,          & ! intent(in)    : mask for the fluxes used in this given state subset
                       ! input/output: data structures
                       model_decisions,   & ! intent(in)    : model decisions
                       type_data,         & ! intent(in)    : type of vegetation and soil
                       attr_data,         & ! intent(in)    : spatial attributes
                       forc_data,         & ! intent(in)    : model forcing data
                       mpar_data,         & ! intent(in)    : model parameters
                       indx_data,         & ! intent(inout) : index data
                       prog_data,         & ! intent(inout) : model prognostic variables for a local HRU
                       diag_data,         & ! intent(inout) : model diagnostic variables for a local HRU
                       flux_data,         & ! intent(inout) : model fluxes for a local HRU
                       deriv_data,        & ! intent(inout) : derivatives in model fluxes w.r.t. relevant state variables
                       bvar_data,         & ! intent(in)    : model variables for the local basin
                       ! output: model control
                       ixSaturation,      & ! intent(inout) : index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       dtMultiplier,      & ! intent(out)   : substep multiplier (-)
                       nSubsteps,         & ! intent(out)   : number of substeps taken for a given split
                       failedMinimumStep, & ! intent(out)   : flag to denote success of substepping for a given split
                       reduceCoupledStep, & ! intent(out)   : flag to denote need to reduce the length of the coupled step
                       tooMuchMelt,       & ! intent(out)   : flag to denote that ice is insufficient to support melt
                       err,message)         ! intent(out)   : error code and error message
 ! ---------------------------------------------------------------------------------------
 ! structure allocations
 USE globalData,only:flux_meta                        ! metadata on the model fluxes
 USE allocspace_module,only:allocLocal                ! allocate local data structures
 ! simulation of fluxes and residuals given a trial state vector
 USE systemSolv_module,only:systemSolv                ! solve the system of equations for one time step
 USE getVectorz_module,only:popStateVec               ! populate the state vector
 USE getVectorz_module,only:varExtract                ! extract variables from the state vector
 USE updateVars_module,only:updateVars                ! update prognostic variables
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                            ! time step (seconds)
 real(dp),intent(in)             :: dt_min                        ! minimum time step (seconds)
 real(dp),intent(in)             :: errTol                        ! error tolerance in the explicit solution
 integer(i4b),intent(in)         :: nState                        ! total number of state variables
 logical(lgt),intent(in)         :: doAdjustTemp                  ! flag to indicate if we adjust the temperature
 logical(lgt),intent(in)         :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall                 ! flag to define the first flux call
 logical(lgt),intent(in)         :: explicitEuler                 ! flag to denote computing the explicit Euler solution
 logical(lgt),intent(in)         :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 logical(lgt),intent(in)         :: fluxMask(:)                   ! flags to denote if the flux is calculated in the given state subset
 ! input/output: data structures
 type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
 type(var_i),intent(in)          :: type_data                     ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                     ! spatial attributes
 type(var_d),intent(in)          :: forc_data                     ! model forcing data
 type(var_dlength),intent(in)    :: mpar_data                     ! model parameters
 type(var_ilength),intent(inout) :: indx_data                     ! indices for a local HRU
 type(var_dlength),intent(inout) :: prog_data                     ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                     ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                     ! model fluxes for a local HRU
 type(var_dlength),intent(inout) :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
 type(var_dlength),intent(in)    :: bvar_data                     ! model variables for the local basin
 ! output: model control
 integer(i4b),intent(inout)      :: ixSaturation                  ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(dp),intent(out)            :: dtMultiplier                  ! substep multiplier (-)
 integer(i4b),intent(out)        :: nSubsteps                     ! number of substeps taken for a given split
 logical(lgt),intent(out)        :: failedMinimumStep             ! flag to denote success of substepping for a given split
 logical(lgt),intent(out)        :: reduceCoupledStep             ! flag to denote need to reduce the length of the coupled step
 logical(lgt),intent(out)        :: tooMuchMelt                   ! flag to denote that ice is insufficient to support melt
 integer(i4b),intent(out)        :: err                           ! error code
 character(*),intent(out)        :: message                       ! error message
 ! ---------------------------------------------------------------------------------------
 ! * general local variables
 ! ---------------------------------------------------------------------------------------
 ! error control
 character(LEN=256)              :: cmessage                      ! error message of downwind routine
 ! time stepping
 real(dp)                        :: dtSum                         ! sum of time from successful steps (seconds)
 real(dp)                        :: dt_wght                       ! weight given to a given flux calculation
 real(dp)                        :: dtSubstep                     ! length of a substep (s) 
 integer(i4b)                    :: iVar                          ! index of variables in data structures
 ! adaptive sub-stepping for the explicit solution
 logical(lgt)                    :: failedSubstep                 ! flag to denote success of substepping for a given split
 real(dp)                        :: explicitError                 ! error in the explicit solution
 real(dp),parameter              :: safety=0.85_dp                ! safety factor in adaptive sub-stepping
 real(dp),parameter              :: reduceMin=0.1_dp              ! mimimum factor that time step is reduced
 real(dp),parameter              :: increaseMax=4.0_dp            ! maximum factor that time step is increased
 ! adaptive sub-stepping for the implicit solution
 integer(i4b)                    :: niter                         ! number of iterations taken
 integer(i4b),parameter          :: n_inc=5                       ! minimum number of iterations to increase time step
 integer(i4b),parameter          :: n_dec=15                      ! maximum number of iterations to decrease time step
 real(dp),parameter              :: F_inc = 1.25_dp               ! factor used to increase time step
 real(dp),parameter              :: F_dec = 0.90_dp               ! factor used to decrease time step
 ! state and flux vectors
 real(dp)                        :: untappedMelt(nState)          ! un-tapped melt energy (J m-3 s-1)
 real(dp)                        :: stateVecInit(nState)          ! initial state vector (mixed units)
 real(dp)                        :: stateVecTrial(nState)         ! trial state vector (mixed units)
 type(var_dlength)               :: flux_temp                     ! temporary model fluxes
 ! flags
 logical(lgt)                    :: waterBalanceError              ! flag to denote that there is a water balance error
 logical(lgt)                    :: nrgFluxModified               ! flag to denote that the energy fluxes were modified
 ! energy fluxes
 real(dp)                        :: sumCanopyEvaporation          ! sum of canopy evaporation/condensation (kg m-2 s-1)
 real(dp)                        :: sumLatHeatCanopyEvap          ! sum of latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
 real(dp)                        :: sumSenHeatCanopy              ! sum of sensible heat flux from the canopy to the canopy air space (W m-2)
 real(dp)                        :: sumSoilCompress 
 real(dp),allocatable            :: sumLayerCompress(:) 
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 globalVars: associate(&
 ! number of layers
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):    [i4b]    number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):    [i4b]    number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):    [i4b]    total number of layers
 mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):    [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
 ! model state variables (vegetation canopy)
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(inout): [dp]     temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(inout): [dp]     temperature of the vegetation canopy (K)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(inout): [dp]     mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(inout): [dp]     mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(inout): [dp]     mass of total water on the vegetation canopy (kg m-2)
 ! model state variables (snow and soil domains)
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(inout): [dp(:)]  temperature of each snow/soil layer (K)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of ice (-)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of liquid water (-)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(inout): [dp(:)]  matric head (m)
 mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat        & ! intent(inout): [dp(:)]  matric potential of liquid water (m)
 )  ! end association with variables in the data structures
 ! *********************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************
 ! Procedure starts here

 ! initialize error control
 err=0; message='varSubstep/'

 ! initialize flag for the success of the substepping
 failedMinimumStep=.false.

 ! initialize the length of the substep
 dtSubstep    = dt

 ! allocate space for the temporary model flux structure
 call allocLocal(flux_meta(:),flux_temp,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! initialize the model fluxes (some model fluxes are not computed in the iterations)
 do iVar=1,size(flux_data%var)
  flux_temp%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:)
 end do

 ! initialize the total energy fluxes (modified in updateProg)
 sumCanopyEvaporation = 0._dp  ! canopy evaporation/condensation (kg m-2 s-1)
 sumLatHeatCanopyEvap = 0._dp  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
 sumSenHeatCanopy     = 0._dp  ! sensible heat flux from the canopy to the canopy air space (W m-2)
 sumSoilCompress      = 0._dp  ! total soil compression
 allocate(sumLayerCompress(nSoil)); sumLayerCompress = 0._dp ! soil compression by layer 

 ! initialize subStep
 dtSum     = 0._dp  ! keep track of the portion of the time step that is completed
 nSubsteps = 0

 ! loop through substeps
 ! NOTE: continuous do statement with exit clause
 substeps: do

  ! initialize error control
  err=0; message='varSubstep/'

  ! increment substep
  nSubsteps = nSubsteps+1

  !print*, '** new substep'
  !print*, 'scalarCanopyIce  = ', prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)
  !print*, 'scalarCanopyTemp = ', prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)

  ! -----
  ! * populate state vectors...
  ! ---------------------------

  ! initialize state vectors
  call popStateVec(&
                   ! input
                   nState,                           & ! intent(in):    number of desired state variables
                   prog_data,                        & ! intent(in):    model prognostic variables for a local HRU
                   diag_data,                        & ! intent(in):    model diagnostic variables for a local HRU
                   indx_data,                        & ! intent(in):    indices defining model states and layers
                   ! output
                   stateVecInit,                     & ! intent(out):   initial model state vector (mixed units)
                   err,cmessage)                       ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! -----
  ! * iterative solution...
  ! -----------------------

  ! solve the system of equations for a given state subset
  call systemSolv(&
                  ! input: model control
                  dtSubstep,         & ! intent(in):    time step (s)
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
                  ! output: model control
                  deriv_data,        & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                  ixSaturation,      & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                  untappedMelt,      & ! intent(out):   un-tapped melt energy (J m-3 s-1)
                  stateVecTrial,     & ! intent(out):   updated state vector
                  explicitError,     & ! intent(out):   error in the explicit solution
                  reduceCoupledStep, & ! intent(out):   flag to reduce the length of the coupled step
                  tooMuchMelt,       & ! intent(out):   flag to denote that ice is insufficient to support melt
                  niter,             & ! intent(out):   number of iterations taken
                  err,cmessage)        ! intent(out):   error code and error message
  if(err/=0)then
   message=trim(message)//trim(cmessage)
   if(err>0) return
  endif

  ! if too much melt or need to reduce length of the coupled step then return
  ! NOTE: need to go all the way back to coupled_em and merge snow layers, as all splitting operations need to occur with the same layer geometry
  if(tooMuchMelt .or. reduceCoupledStep) return

  ! identify failure
  failedSubstep = (err<0)
  if(explicitEuler .and. explicitError > errTol) failedSubstep=.true.

  ! check
  if(globalPrintFlag)then
   print*, 'niter, failedSubstep, dtSubstep = ', niter, failedSubstep, dtSubstep
   print*, trim(cmessage)
  endif

  ! reduce step based on failure
  if(failedSubstep)then

   ! get time step reduction
   ! solver failure
   if(err<0)then
    err=0; message='varSubstep/'  ! recover from failed convergence 
    dtMultiplier  = 0.5_dp        ! system failure: step halving
   ! large error in explicit solution
   elseif(explicitEuler .and. explicitError > errTol)then
    dtMultiplier  = max(safety*sqrt(errTol/explicitError), reduceMin)  ! intolerable errors: MUST be explicit Euler
   ! nothing else defined
   else
    message=trim(message)//'unknown failure'
    err=20; return
   endif

  ! successful step
  else

   ! ** implicit Euler: adjust step length based on iteration count
   if(.not.explicitEuler)then
    if(niter<n_inc)then
     dtMultiplier = F_inc
    elseif(niter>n_dec)then
     dtMultiplier = F_dec
    else
     dtMultiplier = 1._dp
    endif

   ! ** explicit Euler: adjust step based on error estimate
   else
    dtMultiplier  = min(safety*sqrt(errTol/explicitError), increaseMax)
   endif  ! switch between explicit and implicit Euler
  endif  ! switch between failure and success

  ! check if we failed the substep
  if(failedSubstep)then
  
   ! check that the substep is greater than the minimum step
   if(dtSubstep*dtMultiplier<dt_min)then

    ! --> if not explicit Euler, then exit and try another solution method
    if(.not.explicitEuler)then
     failedMinimumStep=.true.
     exit subSteps
    endif

   else ! step is still OK
    dtSubstep = dtSubstep*dtMultiplier  
    cycle subSteps
   endif  ! if step is less than the minimum

  endif  ! if failed the substep

  ! -----
  ! * update model fluxes...
  ! ------------------------

  ! NOTE: if we get to here then we are accepting the step

  ! NOTE: we get to here if iterations are successful
  if(err/=0)then
   message=trim(message)//'expect err=0 if updating fluxes'
   return
  endif

  ! update prognostic variables
  call updateProg(dtSubstep,nSnow,nSoil,nLayers,doAdjustTemp,explicitEuler,computeVegFlux,untappedMelt,stateVecTrial, & ! input: model control
                  mpar_data,indx_data,flux_temp,prog_data,diag_data,deriv_data,                                       & ! input-output: data structures
                  waterBalanceError,nrgFluxModified,err,cmessage)                                                       ! output: flags and error control
  if(err/=0)then
   message=trim(message)//trim(cmessage)
   if(err>0) return
  endif

  ! if water balance error then reduce the length of the coupled step
  if(waterBalanceError)then
   reduceCoupledStep=.true.
   return
  endif

  if(globalPrintFlag)&
  print*, trim(cmessage)//': dt = ', dtSubstep

  ! recover from errors in prognostic update
  if(err<0)then

   ! modify step
   err=0  ! error recovery
   dtSubstep = dtSubstep/2._dp

   ! check minimum: fail minimum step if there is an error in the update
   if(dtSubstep<dt_min)then
    failedMinimumStep=.true.
    exit subSteps
   ! minimum OK -- try again
   else
    cycle substeps
   endif

  endif  ! if errors in prognostic update

  ! get the total energy fluxes (modified in updateProg)
  if(nrgFluxModified .or. indx_data%var(iLookINDEX%ixVegNrg)%dat(1)/=integerMissing)then
   sumCanopyEvaporation = sumCanopyEvaporation + dtSubstep*flux_temp%var(iLookFLUX%scalarCanopyEvaporation)%dat(1) ! canopy evaporation/condensation (kg m-2 s-1)
   sumLatHeatCanopyEvap = sumLatHeatCanopyEvap + dtSubstep*flux_temp%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1) ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
   sumSenHeatCanopy     = sumSenHeatCanopy     + dtSubstep*flux_temp%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)     ! sensible heat flux from the canopy to the canopy air space (W m-2)
  else
   sumCanopyEvaporation = sumCanopyEvaporation + dtSubstep*flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1) ! canopy evaporation/condensation (kg m-2 s-1)
   sumLatHeatCanopyEvap = sumLatHeatCanopyEvap + dtSubstep*flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1) ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
   sumSenHeatCanopy     = sumSenHeatCanopy     + dtSubstep*flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)     ! sensible heat flux from the canopy to the canopy air space (W m-2)
  endif  ! if energy fluxes were modified

  ! get the total soil compressions
  if (count(indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat/=integerMissing)>0) then
   sumSoilCompress = sumSoilCompress + diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1) ! total soil compression
   sumLayerCompress = sumLayerCompress + diag_data%var(iLookDIAG%mLayerCompress)%dat ! soil compression in layers
  endif

  ! print progress
  if(globalPrintFlag)&
  write(*,'(a,1x,3(f13.2,1x))') 'updating: dtSubstep, dtSum, dt = ', dtSubstep, dtSum, dt

  ! increment fluxes
  dt_wght = dtSubstep/dt ! (define weight applied to each splitting operation)
  do iVar=1,size(flux_meta)
   if(fluxMask(iVar)) flux_data%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:)*dt_wght
  end do

  ! ------------------------------------------------------
  ! ------------------------------------------------------

  ! increment sub-step
  dtSum = dtSum + dtSubstep
  !print*, 'dtSum, dtSubstep, dt = ', dtSum, dtSubstep, dt

  ! check that we have completed the sub-step
  if(dtSum >= dt-verySmall)then
   failedMinimumStep=.false.
   exit subSteps
  endif

  ! adjust length of the sub-step (make sure that we don't exceed the step)
  dtSubstep = min(dt - dtSum, max(dtSubstep*dtMultiplier, dt_min) )

 end do substeps  ! time steps for variable-dependent sub-stepping

 ! save the energy fluxes
 flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1) = sumCanopyEvaporation /dt      ! canopy evaporation/condensation (kg m-2 s-1)
 flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1) = sumLatHeatCanopyEvap /dt      ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
 flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)     = sumSenHeatCanopy     /dt      ! sensible heat flux from the canopy to the canopy air space (W m-2)

 ! save the soil compression diagnostics 
 diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1) = sumSoilCompress
 diag_data%var(iLookDIAG%mLayerCompress)%dat = sumLayerCompress
 deallocate(sumLayerCompress)

 ! end associate statements
 end associate globalVars

 end subroutine varSubstep


 ! **********************************************************************************************************
 ! private subroutine updateProg: update prognostic variables
 ! **********************************************************************************************************
 subroutine updateProg(dt,nSnow,nSoil,nLayers,doAdjustTemp,explicitEuler,computeVegFlux,untappedMelt,stateVecTrial, & ! input: model control
                       mpar_data,indx_data,flux_data,prog_data,diag_data,deriv_data,                                & ! input-output: data structures
                       waterBalanceError,nrgFluxModified,err,message)                                                 ! output: flags and error control
 USE getVectorz_module,only:varExtract                             ! extract variables from the state vector
 USE updateVars_module,only:updateVars                             ! update prognostic variables
 implicit none
 ! model control
 real(dp)         ,intent(in)    :: dt                             ! time step (s)
 integer(i4b)     ,intent(in)    :: nSnow                          ! number of snow layers
 integer(i4b)     ,intent(in)    :: nSoil                          ! number of soil layers
 integer(i4b)     ,intent(in)    :: nLayers                        ! total number of layers
 logical(lgt)     ,intent(in)    :: doAdjustTemp                   ! flag to indicate if we adjust the temperature
 logical(lgt)     ,intent(in)    :: explicitEuler                  ! flag to denote computing the explicit Euler solution
 logical(lgt)     ,intent(in)    :: computeVegFlux                 ! flag to compute the vegetation flux
 real(dp)         ,intent(in)    :: untappedMelt(:)                ! un-tapped melt energy (J m-3 s-1)
 real(dp)         ,intent(in)    :: stateVecTrial(:)               ! trial state vector (mixed units)
 ! data structures 
 type(var_dlength),intent(in)    :: mpar_data                      ! model parameters
 type(var_ilength),intent(in)    :: indx_data                      ! indices for a local HRU
 type(var_dlength),intent(inout) :: flux_data                      ! model fluxes for a local HRU
 type(var_dlength),intent(inout) :: prog_data                      ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                      ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: deriv_data                     ! derivatives in model fluxes w.r.t. relevant state variables
 ! flags and error control
 logical(lgt)     ,intent(out)   :: waterBalanceError              ! flag to denote that there is a water balance error
 logical(lgt)     ,intent(out)   :: nrgFluxModified                ! flag to denote that the energy fluxes were modified
 integer(i4b)     ,intent(out)   :: err                            ! error code
 character(*)     ,intent(out)   :: message                        ! error message
 ! ==================================================================================================================
 ! general
 integer(i4b)                    :: iState                         ! index of model state variable
 integer(i4b)                    :: ixSubset                       ! index within the state subset
 integer(i4b)                    :: ixFullVector                   ! index within full state vector
 integer(i4b)                    :: ixControlIndex                 ! index within a given domain
 real(dp)                        :: volMelt                        ! volumetric melt (kg m-3)
 real(dp),parameter              :: verySmall=epsilon(1._dp)*2._dp ! a very small number (deal with precision issues)
 ! mass balance
 logical(lgt),parameter          :: checkMassBalance=.true.        ! flag to check the mass balance
 real(dp)                        :: canopyBalance0,canopyBalance1  ! canopy storage at start/end of time step
 real(dp)                        :: soilBalance0,soilBalance1      ! soil storage at start/end of time step
 real(dp)                        :: vertFlux                       ! change in storage due to vertical fluxes
 real(dp)                        :: tranSink,baseSink,compSink     ! change in storage due to sink terms
 real(dp)                        :: liqError                       ! water balance error
 real(dp)                        :: fluxNet                        ! net water fluxes (kg m-2 s-1)
 real(dp)                        :: superflousWat                  ! superflous water used for evaporation (kg m-2 s-1)
 real(dp)                        :: superflousNrg                  ! superflous energy that cannot be used for evaporation (W m-2 [J m-2 s-1])
 character(LEN=256)              :: cmessage                       ! error message of downwind routine
 ! trial state variables
 real(dp)                        :: scalarCanairTempTrial          ! trial value for temperature of the canopy air space (K)
 real(dp)                        :: scalarCanopyTempTrial          ! trial value for temperature of the vegetation canopy (K)
 real(dp)                        :: scalarCanopyWatTrial           ! trial value for liquid water storage in the canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerTempTrial                ! trial vector for temperature of layers in the snow and soil domains (K)
 real(dp),dimension(nLayers)     :: mLayerVolFracWatTrial          ! trial vector for volumetric fraction of total water (-)
 real(dp),dimension(nSoil)       :: mLayerMatricHeadTrial          ! trial vector for total water matric potential (m)
 real(dp),dimension(nSoil)       :: mLayerMatricHeadLiqTrial       ! trial vector for liquid water matric potential (m)
 ! diagnostic variables
 real(dp)                        :: scalarCanopyLiqTrial           ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 real(dp)                        :: scalarCanopyIceTrial           ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerVolFracLiqTrial          ! trial vector for volumetric fraction of liquid water (-)
 real(dp),dimension(nLayers)     :: mLayerVolFracIceTrial          ! trial vector for volumetric fraction of ice (-)
 ! -------------------------------------------------------------------------------------------------------------------

 ! -------------------------------------------------------------------------------------------------------------------
 ! point to flux variables in the data structure
 associate(&
 ! get indices for mass balance
 ixVegHyd                  => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)                  ,& ! intent(in)    : [i4b]    index of canopy hydrology state variable (mass)
 ixSoilOnlyHyd             => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat                ,& ! intent(in)    : [i4b(:)] index in the state subset for hydrology state variables in the soil domain
 ! get indices for the un-tapped melt
 ixNrgOnly                 => indx_data%var(iLookINDEX%ixNrgOnly)%dat                    ,& ! intent(in)    : [i4b(:)] list of indices for all energy states
 ixDomainType              => indx_data%var(iLookINDEX%ixDomainType)%dat                 ,& ! intent(in)    : [i4b(:)] indices defining the domain of the state (iname_veg, iname_snow, iname_soil)
 ixControlVolume           => indx_data%var(iLookINDEX%ixControlVolume)%dat              ,& ! intent(in)    : [i4b(:)] index of the control volume for different domains (veg, snow, soil)
 ixMapSubset2Full          => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat             ,& ! intent(in)    : [i4b(:)] [state subset] list of indices of the full state vector in the state subset
 ! water fluxes
 scalarRainfall            => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)             ,& ! intent(in)    : [dp]     rainfall rate (kg m-2 s-1)
 scalarThroughfallRain     => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)      ,& ! intent(in)    : [dp]     rain reaches ground without touching the canopy (kg m-2 s-1)
 scalarCanopyEvaporation   => flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)    ,& ! intent(in)    : [dp]     canopy evaporation/condensation (kg m-2 s-1)
 scalarCanopyTranspiration => flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1)  ,& ! intent(in)    : [dp]     canopy transpiration (kg m-2 s-1)
 scalarCanopyLiqDrainage   => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)    ,& ! intent(in)    : [dp]     drainage liquid water from vegetation canopy (kg m-2 s-1)
 iLayerLiqFluxSoil         => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat             ,& ! intent(in)    : [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
 mLayerTranspire           => flux_data%var(iLookFLUX%mLayerTranspire)%dat               ,& ! intent(in)    : [dp(:)]  transpiration loss from each soil layer (m s-1)
 mLayerBaseflow            => flux_data%var(iLookFLUX%mLayerBaseflow)%dat                ,& ! intent(in)    : [dp(:)]  baseflow from each soil layer (m s-1)
 mLayerCompress            => diag_data%var(iLookDIAG%mLayerCompress)%dat                ,& ! intent(in)    : [dp(:)]  change in storage associated with compression of the soil matrix (-)
 scalarCanopySublimation   => flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1)    ,& ! intent(in)    : [dp]     sublimation of ice from the vegetation canopy (kg m-2 s-1)
 scalarSnowSublimation     => flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1)      ,& ! intent(in)    : [dp]     sublimation of ice from the snow surface (kg m-2 s-1)
 ! energy fluxes
 scalarLatHeatCanopyEvap   => flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)    ,& ! intent(in)    : [dp]     latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
 scalarSenHeatCanopy       => flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)        ,& ! intent(in)    : [dp]     sensible heat flux from the canopy to the canopy air space (W m-2)
 ! domain depth
 canopyDepth               => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)          ,& ! intent(in)    : [dp   ]  canopy depth (m)
 mLayerDepth               => prog_data%var(iLookPROG%mLayerDepth)%dat                   ,& ! intent(in)    : [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
 ! model state variables (vegetation canopy)
 scalarCanairTemp          => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)           ,& ! intent(inout) : [dp]     temperature of the canopy air space (K)
 scalarCanopyTemp          => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)           ,& ! intent(inout) : [dp]     temperature of the vegetation canopy (K)
 scalarCanopyIce           => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)            ,& ! intent(inout) : [dp]     mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiq           => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)            ,& ! intent(inout) : [dp]     mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyWat           => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)            ,& ! intent(inout) : [dp]     mass of total water on the vegetation canopy (kg m-2)
 ! model state variables (snow and soil domains)
 mLayerTemp                => prog_data%var(iLookPROG%mLayerTemp)%dat                    ,& ! intent(inout) : [dp(:)]  temperature of each snow/soil layer (K)
 mLayerVolFracIce          => prog_data%var(iLookPROG%mLayerVolFracIce)%dat              ,& ! intent(inout) : [dp(:)]  volumetric fraction of ice (-)
 mLayerVolFracLiq          => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat              ,& ! intent(inout) : [dp(:)]  volumetric fraction of liquid water (-)
 mLayerVolFracWat          => prog_data%var(iLookPROG%mLayerVolFracWat)%dat              ,& ! intent(inout) : [dp(:)]  volumetric fraction of total water (-)
 mLayerMatricHead          => prog_data%var(iLookPROG%mLayerMatricHead)%dat              ,& ! intent(inout) : [dp(:)]  matric head (m)
 mLayerMatricHeadLiq       => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat           ,& ! intent(inout) : [dp(:)]  matric potential of liquid water (m)
 ! error tolerance
 absConvTol_liquid         => mpar_data%var(iLookPARAM%absConvTol_liquid)%dat(1)          & ! intent(in)    : [dp]     absolute convergence tolerance for vol frac liq water (-)
 ) ! associating flux variables in the data structure
 ! -------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='updateProg/'

 ! initialize water balance error
 waterBalanceError=.false.

 ! get storage at the start of the step
 canopyBalance0 = merge(scalarCanopyWat, realMissing, computeVegFlux)
 soilBalance0   = sum( (mLayerVolFracLiq(nSnow+1:nLayers)  + mLayerVolFracIce(nSnow+1:nLayers)  )*mLayerDepth(nSnow+1:nLayers) )

 ! -----
 ! * update states...
 ! ------------------

 ! extract states from the state vector
 call varExtract(&
                 ! input
                 stateVecTrial,            & ! intent(in):    model state vector (mixed units)
                 diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                 prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                 indx_data,                & ! intent(in):    indices defining model states and layers
                 ! output: variables for the vegetation canopy
                 scalarCanairTempTrial,    & ! intent(out):   trial value of canopy air temperature (K)
                 scalarCanopyTempTrial,    & ! intent(out):   trial value of canopy temperature (K)
                 scalarCanopyWatTrial,     & ! intent(out):   trial value of canopy total water (kg m-2)
                 scalarCanopyLiqTrial,     & ! intent(out):   trial value of canopy liquid water (kg m-2)
                 scalarCanopyIceTrial,     & ! intent(out):   trial value of canopy ice content (kg m-2)
                 ! output: variables for the snow-soil domain
                 mLayerTempTrial,          & ! intent(out):   trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,    & ! intent(out):   trial vector of volumetric total water content (-)
                 mLayerVolFracLiqTrial,    & ! intent(out):   trial vector of volumetric liquid water content (-)
                 mLayerVolFracIceTrial,    & ! intent(out):   trial vector of volumetric ice water content (-)
                 mLayerMatricHeadTrial,    & ! intent(out):   trial vector of total water matric potential (m)
                 mLayerMatricHeadLiqTrial, & ! intent(out):   trial vector of liquid water matric potential (m)
                 ! output: error control
                 err,cmessage)               ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

 !print*, 'after varExtract: scalarCanopyTempTrial =', scalarCanopyTempTrial   ! trial value of canopy temperature (K)
 !print*, 'after varExtract: scalarCanopyWatTrial  =', scalarCanopyWatTrial    ! trial value of canopy total water (kg m-2)
 !print*, 'after varExtract: scalarCanopyLiqTrial  =', scalarCanopyLiqTrial    ! trial value of canopy liquid water (kg m-2)
 !print*, 'after varExtract: scalarCanopyIceTrial  =', scalarCanopyIceTrial    ! trial value of canopy ice content (kg m-2)

 ! update diagnostic variables
 call updateVars(&
                 ! input
                 doAdjustTemp,             & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                 explicitEuler,            & ! intent(in):    flag to denote computing the explicit Euler solution
                 mpar_data,                & ! intent(in):    model parameters for a local HRU
                 indx_data,                & ! intent(in):    indices defining model states and layers
                 prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                 diag_data,                & ! intent(inout): model diagnostic variables for a local HRU
                 deriv_data,               & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! output: variables for the vegetation canopy
                 scalarCanopyTempTrial,    & ! intent(inout): trial value of canopy temperature (K)
                 scalarCanopyWatTrial,     & ! intent(inout): trial value of canopy total water (kg m-2)
                 scalarCanopyLiqTrial,     & ! intent(inout): trial value of canopy liquid water (kg m-2)
                 scalarCanopyIceTrial,     & ! intent(inout): trial value of canopy ice content (kg m-2)
                 ! output: variables for the snow-soil domain
                 mLayerTempTrial,          & ! intent(inout): trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,    & ! intent(inout): trial vector of volumetric total water content (-)
                 mLayerVolFracLiqTrial,    & ! intent(inout): trial vector of volumetric liquid water content (-)
                 mLayerVolFracIceTrial,    & ! intent(inout): trial vector of volumetric ice water content (-)
                 mLayerMatricHeadTrial,    & ! intent(inout): trial vector of total water matric potential (m)
                 mLayerMatricHeadLiqTrial, & ! intent(inout): trial vector of liquid water matric potential (m)
                 ! output: error control
                 err,cmessage)               ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

 !print*, 'after updateVars: scalarCanopyTempTrial =', scalarCanopyTempTrial   ! trial value of canopy temperature (K)
 !print*, 'after updateVars: scalarCanopyWatTrial  =', scalarCanopyWatTrial    ! trial value of canopy total water (kg m-2)
 !print*, 'after updateVars: scalarCanopyLiqTrial  =', scalarCanopyLiqTrial    ! trial value of canopy liquid water (kg m-2)
 !print*, 'after updateVars: scalarCanopyIceTrial  =', scalarCanopyIceTrial    ! trial value of canopy ice content (kg m-2)

 ! -----
 ! * check mass balance...
 ! -----------------------

 ! NOTE: should not need to do this, since mass balance is checked in the solver
 if(checkMassBalance)then

  ! check mass balance for the canopy
  if(ixVegHyd/=integerMissing)then

   ! handle cases where fluxes empty the canopy
   fluxNet = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
   if(-fluxNet*dt > canopyBalance0)then

    ! --> first add water
    canopyBalance1 = canopyBalance0 + (scalarRainfall - scalarThroughfallRain)*dt

    ! --> next, remove canopy evaporation -- put the unsatisfied evap into sensible heat
    canopyBalance1 = canopyBalance1 + scalarCanopyEvaporation*dt
    if(canopyBalance1 < 0._dp)then
     ! * get superfluous water and energy
     superflousWat = -canopyBalance1/dt     ! kg m-2 s-1
     superflousNrg = superflousWat*LH_vap   ! W m-2 (J m-2 s-1)
     ! * update fluxes and states
     canopyBalance1          = 0._dp
     scalarCanopyEvaporation = scalarCanopyEvaporation + superflousWat
     scalarLatHeatCanopyEvap = scalarLatHeatCanopyEvap + superflousNrg
     scalarSenHeatCanopy     = scalarSenHeatCanopy - superflousNrg
    endif

    ! --> next, remove canopy drainage
    canopyBalance1 = canopyBalance1 - scalarCanopyLiqDrainage*dt
    if(canopyBalance1 < 0._dp)then
     superflousWat            = -canopyBalance1/dt     ! kg m-2 s-1
     canopyBalance1          = 0._dp
     scalarCanopyLiqDrainage = scalarCanopyLiqDrainage + superflousWat 
    endif

    ! update the trial state
    scalarCanopyWatTrial = canopyBalance1
 
    ! set the modification flag
    nrgFluxModified = .true.

   else
    canopyBalance1  = canopyBalance0 + fluxNet*dt
    nrgFluxModified = .false.
   endif  ! cases where fluxes empty the canopy

   ! check the mass balance
   fluxNet  = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
   liqError = (canopyBalance0 + fluxNet*dt) - scalarCanopyWatTrial
   !write(*,'(a,1x,f20.10)') 'dt = ', dt
   !write(*,'(a,1x,f20.10)') 'scalarCanopyWatTrial       = ', scalarCanopyWatTrial
   !write(*,'(a,1x,f20.10)') 'canopyBalance0             = ', canopyBalance0
   !write(*,'(a,1x,f20.10)') 'canopyBalance1             = ', canopyBalance1
   !write(*,'(a,1x,f20.10)') 'scalarRainfall*dt          = ', scalarRainfall*dt
   !write(*,'(a,1x,f20.10)') 'scalarCanopyLiqDrainage*dt = ', scalarCanopyLiqDrainage*dt
   !write(*,'(a,1x,f20.10)') 'scalarCanopyEvaporation*dt = ', scalarCanopyEvaporation*dt
   !write(*,'(a,1x,f20.10)') 'scalarThroughfallRain*dt   = ', scalarThroughfallRain*dt
   !write(*,'(a,1x,f20.10)') 'liqError                   = ', liqError 
   if(abs(liqError) > absConvTol_liquid*10._dp)then  ! *10 because of precision issues
    waterBalanceError = .true.
    return
   endif  ! if there is a water balance error
  endif  ! if veg canopy
  
  ! check mass balance for soil
  ! NOTE: fatal errors, though possible to recover using negative error codes
  if(count(ixSoilOnlyHyd/=integerMissing)>0)then
   soilBalance1 = sum( (mLayerVolFracLiqTrial(nSnow+1:nLayers) + mLayerVolFracIceTrial(nSnow+1:nLayers) )*mLayerDepth(nSnow+1:nLayers) )
   vertFlux     = -(iLayerLiqFluxSoil(nSoil) - iLayerLiqFluxSoil(0))*dt  ! m s-1 --> m
   tranSink     = sum(mLayerTranspire)*dt                                ! m s-1 --> m
   baseSink     = sum(mLayerBaseflow)*dt                                 ! m s-1 --> m
   compSink     = sum(mLayerCompress(1:nSoil) * mLayerDepth(nSnow+1:nLayers) ) ! dimensionless --> m
   liqError     = soilBalance1 - (soilBalance0 + vertFlux + tranSink - baseSink - compSink)
   !write(*,'(a,1x,f20.10)') 'dt = ', dt
   !write(*,'(a,1x,f20.10)') 'soilBalance0      = ', soilBalance0
   !write(*,'(a,1x,f20.10)') 'soilBalance1      = ', soilBalance1
   !write(*,'(a,1x,f20.10)') 'vertFlux          = ', vertFlux
   !write(*,'(a,1x,f20.10)') 'tranSink          = ', tranSink
   !write(*,'(a,1x,f20.10)') 'baseSink          = ', baseSink
   !write(*,'(a,1x,f20.10)') 'compSink          = ', compSink
   !write(*,'(a,1x,f20.10)') 'liqError          = ', liqError
   !write(*,'(a,1x,f20.10)') 'absConvTol_liquid = ', absConvTol_liquid
   if(abs(liqError) > absConvTol_liquid*10._dp)then   ! *10 because of precision issues
    waterBalanceError = .true.
    return
   endif  ! if there is a water balance error
  endif  ! if hydrology states exist in the soil domain

 endif  ! if checking the mass balance

 ! -----
 ! * remove untapped melt energy (in the explicit Euler method)...
 ! ---------------------------------------------------------------

 ! only work with energy state variables
 if(size(ixNrgOnly)>0)then  ! energy state variables exist

  ! loop through energy state variables
  do iState=1,size(ixNrgOnly)

   ! get index of the control volume within the domain
   ixSubset       = ixNrgOnly(iState)             ! index within the state subset
   ixFullVector   = ixMapSubset2Full(ixSubset)    ! index within full state vector
   ixControlIndex = ixControlVolume(ixFullVector) ! index within a given domain

   ! compute volumetric melt (kg m-3)
   volMelt = dt*untappedMelt(ixSubset)/LH_fus  ! (kg m-3)

   ! update ice content
   select case( ixDomainType(ixFullVector) )
    case(iname_cas);  cycle ! do nothing, since there is no snow stored in the canopy air space
    case(iname_veg);  scalarCanopyIceTrial                        = scalarCanopyIceTrial                        - volMelt*canopyDepth  ! (kg m-2)
    case(iname_snow); mLayerVolFracIceTrial(ixControlIndex)       = mLayerVolFracIceTrial(ixControlIndex)       - volMelt/iden_ice     ! (-)
    case(iname_soil); mLayerVolFracIceTrial(ixControlIndex+nSnow) = mLayerVolFracIceTrial(ixControlIndex+nSnow) - volMelt/iden_water   ! (-)
    case default; err=20; message=trim(message)//'unable to identify domain type [remove untapped melt energy]'; return
   end select

   ! update liquid water content
   select case( ixDomainType(ixFullVector) )
    case(iname_cas);  cycle ! do nothing, since there is no snow stored in the canopy air space
    case(iname_veg);  scalarCanopyLiqTrial                        = scalarCanopyLiqTrial                        + volMelt*canopyDepth  ! (kg m-2)
    case(iname_snow); mLayerVolFracLiqTrial(ixControlIndex)       = mLayerVolFracLiqTrial(ixControlIndex)       + volMelt/iden_water   ! (-)
    case(iname_soil); mLayerVolFracLiqTrial(ixControlIndex+nSnow) = mLayerVolFracLiqTrial(ixControlIndex+nSnow) + volMelt/iden_water   ! (-)
    case default; err=20; message=trim(message)//'unable to identify domain type [remove untapped melt energy]'; return
   end select

  end do  ! looping through energy variables

  ! ========================================================================================================

  ! *** ice

  ! --> check if we removed too much water
  if(scalarCanopyIceTrial < 0._dp  .or. any(mLayerVolFracIceTrial < 0._dp) )then

   ! **
   ! canopy within numerical precision
   if(scalarCanopyIceTrial < 0._dp)then

    if(scalarCanopyIceTrial > -verySmall)then
     scalarCanopyLiqTrial = scalarCanopyLiqTrial - scalarCanopyIceTrial
     scalarCanopyIceTrial = 0._dp
  
    ! encountered an inconsistency: spit the dummy
    else
     print*, 'dt = ', dt
     print*, 'untappedMelt          = ', untappedMelt
     print*, 'untappedMelt*dt       = ', untappedMelt*dt
     print*, 'scalarCanopyiceTrial  = ', scalarCanopyIceTrial
     message=trim(message)//'melted more than the available water'
     err=20; return
    endif  ! (inconsistency)

   endif  ! if checking the canopy

   ! **
   ! snow+soil within numerical precision
   do iState=1,size(mLayerVolFracIceTrial) 

    ! snow layer within numerical precision
    if(mLayerVolFracIceTrial(iState) < 0._dp)then

     if(mLayerVolFracIceTrial(iState) > -verySmall)then
      mLayerVolFracLiqTrial(iState) = mLayerVolFracLiqTrial(iState) - mLayerVolFracIceTrial(iState)
      mLayerVolFracIceTrial(iState) = 0._dp
  
     ! encountered an inconsistency: spit the dummy
     else
      print*, 'dt = ', dt
      print*, 'untappedMelt          = ', untappedMelt
      print*, 'untappedMelt*dt       = ', untappedMelt*dt
      print*, 'mLayerVolFracIceTrial = ', mLayerVolFracIceTrial 
      message=trim(message)//'melted more than the available water'
      err=20; return
     endif  ! (inconsistency)

    endif  ! if checking a snow layer

   end do ! (looping through state variables)

  endif  ! (if we removed too much water)

  ! ========================================================================================================

  ! *** liquid water

  ! --> check if we removed too much water
  if(scalarCanopyLiqTrial < 0._dp  .or. any(mLayerVolFracLiqTrial < 0._dp) )then

   ! **
   ! canopy within numerical precision
   if(scalarCanopyLiqTrial < 0._dp)then

    if(scalarCanopyLiqTrial > -verySmall)then
     scalarCanopyIceTrial = scalarCanopyIceTrial - scalarCanopyLiqTrial
     scalarCanopyLiqTrial = 0._dp

    ! encountered an inconsistency: spit the dummy
    else
     print*, 'dt = ', dt
     print*, 'untappedMelt          = ', untappedMelt
     print*, 'untappedMelt*dt       = ', untappedMelt*dt
     print*, 'scalarCanopyLiqTrial  = ', scalarCanopyLiqTrial
     message=trim(message)//'frozen more than the available water'
     err=20; return
    endif  ! (inconsistency)

   endif  ! checking the canopy
 
   ! **
   ! snow+soil within numerical precision
   do iState=1,size(mLayerVolFracLiqTrial)

    ! snow layer within numerical precision
    if(mLayerVolFracLiqTrial(iState) < 0._dp)then

     if(mLayerVolFracLiqTrial(iState) > -verySmall)then
      mLayerVolFracIceTrial(iState) = mLayerVolFracIceTrial(iState) - mLayerVolFracLiqTrial(iState)
      mLayerVolFracLiqTrial(iState) = 0._dp
  
     ! encountered an inconsistency: spit the dummy
     else
      print*, 'dt = ', dt
      print*, 'untappedMelt          = ', untappedMelt
      print*, 'untappedMelt*dt       = ', untappedMelt*dt
      print*, 'mLayerVolFracLiqTrial = ', mLayerVolFracLiqTrial
      message=trim(message)//'frozen more than the available water'
      err=20; return
     endif  ! (inconsistency)

    endif  ! checking a snow layer

   end do ! (looping through state variables)

  endif  ! (if we removed too much water)

 endif  ! (if energy state variables exist)

 ! -----
 ! * update prognostic variables... 
 ! --------------------------------

 ! build elements of the state vector for the vegetation canopy
 scalarCanairTemp    = scalarCanairTempTrial    ! trial value of canopy air temperature (K)
 scalarCanopyTemp    = scalarCanopyTempTrial    ! trial value of canopy temperature (K)
 scalarCanopyWat     = scalarCanopyWatTrial     ! trial value of canopy total water (kg m-2) 
 scalarCanopyLiq     = scalarCanopyLiqTrial     ! trial value of canopy liquid water (kg m-2)
 scalarCanopyIce     = scalarCanopyIceTrial     ! trial value of canopy ice content (kg m-2) 

 ! build elements of the state vector for the snow+soil domain
 mLayerTemp          = mLayerTempTrial          ! trial vector of layer temperature (K)
 mLayerVolFracWat    = mLayerVolFracWatTrial    ! trial vector of volumetric total water content (-)
 mLayerVolFracLiq    = mLayerVolFracLiqTrial    ! trial vector of volumetric liquid water content (-)
 mLayerVolFracIce    = mLayerVolFracIceTrial    ! trial vector of volumetric ice water content (-)
 mLayerMatricHead    = mLayerMatricHeadTrial    ! trial vector of matric head (m)
 mLayerMatricHeadLiq = mLayerMatricHeadLiqTrial ! trial vector of matric head (m)

 ! end associations to info in the data structures
 end associate

 end subroutine updateProg

end module varSubstep_module
