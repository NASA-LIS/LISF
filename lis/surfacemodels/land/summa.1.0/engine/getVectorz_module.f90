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

module getVectorz_module

! data types
USE nrtype

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! access the global print flag
USE globalData,only:globalPrintFlag

! domain types
USE globalData,only:iname_cas       ! named variables for canopy air space
USE globalData,only:iname_veg       ! named variables for vegetation canopy
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! metadata for information in the data structures
USE globalData,only:indx_meta       ! metadata for the variables in the index structure

! constants
USE multiconst,only:&
                    gravity,      & ! acceleration of gravity              (m s-2)
                    Tfreeze,      & ! temperature at freezing              (K)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (dp)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookDIAG             ! named variables for structure elements
USE var_lookup,only:iLookPROG             ! named variables for structure elements
USE var_lookup,only:iLookDERIV            ! named variables for structure elements
USE var_lookup,only:iLookPARAM            ! named variables for structure elements
USE var_lookup,only:iLookINDEX            ! named variables for structure elements

! provide access to routines to update states
USE updatState_module,only:updateSnow     ! update snow states
USE updatState_module,only:updateSoil     ! update soil states

! provide access to functions for the constitutive functions and derivatives
USE snow_utils_module,only:fracliquid     ! compute the fraction of liquid water (snow) 
USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
USE soil_utils_module,only:dTheta_dTk     ! differentiate the freezing curve w.r.t. temperature (soil)
USE soil_utils_module,only:dTheta_dPsi    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:dPsi_dTheta    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric water content
USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
USE soil_utils_module,only:liquidHead     ! compute the liquid water matric potential 

implicit none
private
public::popStateVec
public::getScaling
public::varExtract

! common variables
real(dp),parameter :: valueMissing=-9999._dp ! missing value

contains

 ! **********************************************************************************************************
 ! public subroutine popStateVec: populate model state vectors 
 ! **********************************************************************************************************
 subroutine popStateVec(&
                        ! input: data structures
                        nState,                  & ! intent(in):    number of desired state variables
                        prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                        indx_data,               & ! intent(in):    indices defining model states and layers
                        ! output
                        stateVec,                & ! intent(out):   model state vector
                        err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: data structures
 integer(i4b),intent(in)         :: nState                 ! number of desired state variables
 type(var_dlength),intent(in)    :: prog_data              ! prognostic variables for a local HRU
 type(var_dlength),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
 type(var_ilength),intent(in)    :: indx_data              ! indices defining model states and layers
 ! output
 real(dp),intent(out)            :: stateVec(:)            ! model state vector (mixed units)
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! state subsets
 integer(i4b)                    :: iState                 ! index of state within the snow+soil domain
 integer(i4b)                    :: iLayer                 ! index of layer within the snow+soil domain
 integer(i4b)                    :: ixStateSubset          ! index within the state subset
 logical(lgt),dimension(nState)  :: stateFlag              ! flag to denote that the state is populated
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 fixedLength: associate(&
 ! model states for the vegetation canopy
 scalarCanairTemp    => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in) : [dp]     temperature of the canopy air space (K)
 scalarCanopyTemp    => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in) : [dp]     temperature of the vegetation canopy (K)
 scalarCanopyWat     => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in) : [dp]     mass of total water on the vegetation canopy (kg m-2)
 scalarCanopyLiq     => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(in) : [dp]     mass of liquid water on the vegetation canopy (kg m-2)
 ! model state variable vectors for the snow-soil layers
 mLayerTemp          => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in) : [dp(:)]  temperature of each snow/soil layer (K)
 mLayerVolFracWat    => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in) : [dp(:)]  volumetric fraction of total water (-)
 mLayerVolFracLiq    => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in) : [dp(:)]  volumetric fraction of liquid water (-)
 mLayerMatricHead    => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in) : [dp(:)]  matric head (m)
 mLayerMatricHeadLiq => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,& ! intent(in) : [dp(:)]  matric potential of liquid water (m)
 ! indices defining specific model states
 ixCasNrg            => indx_data%var(iLookINDEX%ixCasNrg)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy air space energy state variable
 ixVegNrg            => indx_data%var(iLookINDEX%ixVegNrg)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy energy state variable
 ixVegHyd            => indx_data%var(iLookINDEX%ixVegHyd)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy hydrology state variable (mass)
 ! vector of energy and hydrology indices for the snow and soil domains
 ixSnowSoilNrg       => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowSoilHyd       => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 nSnowSoilNrg        => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in) : [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd        => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in) : [i4b]    number of hydrology state variables in the snow+soil domain
 ! type of model state variabless
 ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
 ixHydType           => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in) : [i4b(:)] index of the type of hydrology states in snow+soil domain
 ! number of layers
 nSnow               => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in) : [i4b]    number of snow layers
 nSoil               => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in) : [i4b]    number of soil layers
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in) : [i4b]    total number of layers
 )  ! end association with variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='popStateVec/'

 ! -----
 ! * initialize state vectors...
 ! -----------------------------

 ! initialize flags
 stateFlag(:) = .false.

 ! build the state vector for the temperature of thecanopy air space
 ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer canopy
 do concurrent (iState=1:size(ixCasNrg),ixCasNrg(iState)/=integerMissing)
  stateVec( ixCasNrg(iState) )  = scalarCanairTemp       ! transfer canopy air temperature to the state vector
  stateFlag( ixCasNrg(iState) ) = .true.                 ! flag to denote that the state is populated
 end do

 ! build the state vector for the temperature of the vegetation canopy
 ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer canopy
 do concurrent (iState=1:size(ixVegNrg),ixVegNrg(iState)/=integerMissing)
  stateVec( ixVegNrg(iState) )  = scalarCanopyTemp       ! transfer vegetation temperature to the state vector
  stateFlag( ixVegNrg(iState) ) = .true.                 ! flag to denote that the state is populated
 end do

 ! build the state vector for the water in the vegetation canopy
 ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer canopy
 do concurrent (iState=1:size(ixVegHyd),ixVegHyd(iState)/=integerMissing)
  stateFlag( ixVegHyd(iState) ) = .true.                 ! flag to denote that the state is populated
  select case(ixStateType_subset( ixVegHyd(iState) ))
   case(iname_watCanopy); stateVec( ixVegHyd(iState) )  = scalarCanopyWat        ! transfer total canopy water to the state vector
   case(iname_liqCanopy); stateVec( ixVegHyd(iState) )  = scalarCanopyLiq        ! transfer liquid canopy water to the state vector
   case default; stateFlag( ixVegHyd(iState) ) = .false. ! flag to denote that the state is populated
  end select
 end do

 ! build the energy state vector for the snow and soil domain
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   ixStateSubset            = ixSnowSoilNrg(iLayer)  ! index within the state vector
   stateVec(ixStateSubset)  = mLayerTemp(iLayer)     ! transfer temperature from a layer to the state vector
   stateFlag(ixStateSubset) = .true.                 ! flag to denote that the state is populated
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! build the hydrology state vector for the snow+soil domains
 ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
   ixStateSubset            = ixSnowSoilHyd(iLayer)   ! index within the state vector
   stateFlag(ixStateSubset) = .true.                  ! flag to denote that the state is populated
   select case( ixHydType(iLayer) )
    case(iname_watLayer); stateVec(ixStateSubset) = mLayerVolFracWat(iLayer)           ! total water state variable for snow+soil layers
    case(iname_liqLayer); stateVec(ixStateSubset) = mLayerVolFracLiq(iLayer)           ! liquid water state variable for snow+soil layers
    case(iname_matLayer); stateVec(ixStateSubset) = mLayerMatricHead(iLayer-nSnow)     ! total water matric potential variable for soil layers
    case(iname_lmpLayer); stateVec(ixStateSubset) = mLayerMatricHeadLiq(iLayer-nSnow)  ! liquid matric potential state variable for soil layers
    case default; stateFlag(ixStateSubset) = .false.  ! flag to denote that the state is populated
   end select
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! check that we populated all state variables
 if(count(stateFlag)/=nState)then
  print*, 'stateFlag = ', stateFlag
  message=trim(message)//'some state variables unpopulated'
  err=20; return
 endif

 end associate fixedLength      ! end association to variables in the data structure where vector length does not change
 end subroutine popStateVec


 ! **********************************************************************************************************
 ! public subroutine getScaling: get scale factors 
 ! **********************************************************************************************************
 subroutine getScaling(&
                       ! input: data structures
                       diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                       indx_data,               & ! intent(in):    indices defining model states and layers
                       ! output
                       fScale,                  & ! intent(out):   function scaling vector (mixed units)
                       xScale,                  & ! intent(out):   variable scaling vector (mixed units)
                       sMul,                    & ! intent(out):   multiplier for state vector (used in the residual calculations)
                       dMat,                    & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes) 
                       err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 USE nr_utility_module,only:arth                   ! get a sequence of numbers arth(start, incr, count)
 USE f2008funcs_module,only:findIndex              ! finds the index of the first value within a vector
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: data structures
 type(var_dlength),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
 type(var_ilength),intent(in)    :: indx_data              ! indices defining model states and layers
 ! output: state vectors
 real(dp),intent(out)            :: fScale(:)              ! function scaling vector (mixed units)
 real(dp),intent(out)            :: xScale(:)              ! variable scaling vector (mixed units)
 real(qp),intent(out)            :: sMul(:)    ! NOTE: qp  ! multiplier for state vector (used in the residual calculations)
 real(dp),intent(out)            :: dMat(:)                ! diagonal of the Jacobian matrix (excludes fluxes)
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! scaling parameters
 real(dp),parameter              :: fScaleLiq=0.01_dp      ! func eval: characteristic scale for volumetric liquid water content (-)
 real(dp),parameter              :: fScaleMat=10._dp       ! func eval: characteristic scale for matric head (m)
 real(dp),parameter              :: fScaleNrg=1000000._dp  ! func eval: characteristic scale for energy (J m-3)
 real(dp),parameter              :: xScaleLiq=0.1_dp       ! state var: characteristic scale for volumetric liquid water content (-)
 real(dp),parameter              :: xScaleMat=10._dp       ! state var: characteristic scale for matric head (m)
 real(dp),parameter              :: xScaleTemp=1._dp       ! state var: characteristic scale for temperature (K)
 ! state subsets
 integer(i4b)                    :: iLayer                 ! index of layer within the snow+soil domain
 integer(i4b)                    :: ixStateSubset          ! index within the state subset
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 fixedLength: associate(&
 ! model diagnostic variables
 canopyDepth         => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in):  [dp]     canopy depth (m)
 volHeatCapVeg       => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(in) : [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
 mLayerVolHeatCap    => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,& ! intent(in) : [dp(:)]  bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
 ! indices defining specific model states
 ixCasNrg            => indx_data%var(iLookINDEX%ixCasNrg)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy air space energy state variable
 ixVegNrg            => indx_data%var(iLookINDEX%ixVegNrg)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy energy state variable
 ixVegHyd            => indx_data%var(iLookINDEX%ixVegHyd)%dat                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy hydrology state variable (mass)
 ! vector of energy and hydrology indices for the snow and soil domains
 ixSnowSoilNrg       => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowSoilHyd       => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 nSnowSoilNrg        => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in) : [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd        => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in) : [i4b]    number of hydrology state variables in the snow+soil domain
 ! type of model state variabless
 ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
 ! number of layers
 nSnow               => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in) : [i4b]    number of snow layers
 nSoil               => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in) : [i4b]    number of soil layers
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in) : [i4b]    total number of layers
 )  ! end association with variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='getScaling/'

 ! -----
 ! * define scaling vectors...
 ! ---------------------------

 ! define the function and variable scaling factors for energy
 where(ixStateType_subset==iname_nrgCanair .or. ixStateType_subset==iname_nrgCanopy .or. ixStateType_subset==iname_nrgLayer)
  fScale = 1._dp / fScaleNrg  ! 1/(J m-3)
  xScale = 1._dp  ! K
 endwhere

 ! define the function and variable scaling factors for water on the vegetation canopy
 where(ixStateType_subset==iname_watCanopy .or. ixStateType_subset==iname_liqCanopy)
  fScale = 1._dp / (fScaleLiq*canopyDepth*iden_water)  ! 1/(kg m-2)
  xScale = 1._dp  ! (kg m-2)
 endwhere

 ! define the function and variable scaling factors for water in the snow+soil domain
 where(ixStateType_subset==iname_watLayer .or. ixStateType_subset==iname_liqLayer)
  fScale = 1._dp / fScaleLiq  ! (-)
  xScale = 1._dp  ! (-)
 end where
 
 ! define the function and variable scaling factors for water in the snow+soil domain
 where(ixStateType_subset==iname_matLayer .or. ixStateType_subset==iname_lmpLayer)
  fScale = 1._dp / fScaleLiq  ! (-)
  xScale = 1._dp  ! (m)
 end where

 ! -----
 ! * define components of derivative matrices that are constant over a time step (substep)...
 ! ------------------------------------------------------------------------------------------

 ! define the multiplier for the state vector for residual calculations (vegetation canopy)
 ! NOTE: Use the "where" statement to generalize to multiple canopy layers (currently one canopy layer)
 where(ixStateType_subset==iname_nrgCanair) sMul = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
 where(ixStateType_subset==iname_nrgCanopy) sMul = volHeatCapVeg     ! volumetric heat capacity of the vegetation (J m-3 K-1)
 where(ixStateType_subset==iname_watCanopy) sMul = 1._dp             ! nothing else on the left hand side
 where(ixStateType_subset==iname_liqCanopy) sMul = 1._dp             ! nothing else on the left hand side

 ! compute terms in the Jacobian for vegetation (excluding fluxes)
 ! NOTE: This is computed outside the iteration loop because it does not depend on state variables
 ! NOTE: Energy for vegetation is computed *within* the iteration loop as it includes phase change
 ! NOTE: Use the "where" statement to generalize to multiple canopy layers (currently one canopy layer)
 where(ixStateType_subset==iname_nrgCanair) dMat = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
 where(ixStateType_subset==iname_nrgCanopy) dMat = realMissing       ! populated within the iteration loop 
 where(ixStateType_subset==iname_watCanopy) dMat = 1._dp             ! nothing else on the left hand side
 where(ixStateType_subset==iname_liqCanopy) dMat = 1._dp             ! nothing else on the left hand side

 ! define the energy multiplier and diagonal elements for the state vector for residual calculations (snow-soil domain)
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   ixStateSubset        = ixSnowSoilNrg(iLayer)      ! index within the state vector
   sMul(ixStateSubset)  = mLayerVolHeatCap(iLayer)   ! transfer volumetric heat capacity to the state multiplier
   dMat(ixStateSubset)  = realMissing                ! diagonal element populated within the iteration loop
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! define the hydrology multiplier and diagonal elements for the state vector for residual calculations (snow-soil domain)
 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   ixStateSubset        = ixSnowSoilHyd(iLayer)      ! index within the state vector
   sMul(ixStateSubset)  = 1._dp                      ! state multiplier = 1 (nothing else on the left-hand-side) 
   dMat(ixStateSubset)  = 1._dp                      ! diagonal element = 1 (nothing else on the left-hand-side) 
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! ------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------
 
 end associate fixedLength      ! end association to variables in the data structure where vector length does not change
 end subroutine getScaling



 ! **********************************************************************************************************
 ! public subroutine varExtract: extract variables from the state vector and compute diagnostic variables
 ! **********************************************************************************************************
 subroutine varExtract(&
                       ! input
                       stateVec,                                  & ! intent(in):    model state vector (mixed units)
                       diag_data,                                 & ! intent(in):    model diagnostic variables for a local HRU
                       prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
                       ! output: variables for the vegetation canopy
                       scalarCanairTempTrial,                     & ! intent(out):   trial value of canopy air temperature (K)
                       scalarCanopyTempTrial,                     & ! intent(out):   trial value of canopy temperature (K)
                       scalarCanopyWatTrial,                      & ! intent(out):   trial value of canopy total water (kg m-2)
                       scalarCanopyLiqTrial,                      & ! intent(out):   trial value of canopy liquid water (kg m-2)
                       scalarCanopyIceTrial,                      & ! intent(out):   trial value of canopy ice content (kg m-2)
                       ! output: variables for the snow-soil domain
                       mLayerTempTrial,                           & ! intent(out):   trial vector of layer temperature (K)
                       mLayerVolFracWatTrial,                     & ! intent(out):   trial vector of volumetric total water content (-) 
                       mLayerVolFracLiqTrial,                     & ! intent(out):   trial vector of volumetric liquid water content (-) 
                       mLayerVolFracIceTrial,                     & ! intent(out):   trial vector of volumetric ice water content (-) 
                       mLayerMatricHeadTrial,                     & ! intent(out):   trial vector of total water matric potential (m)
                       mLayerMatricHeadLiqTrial,                  & ! intent(out):   trial vector of liquid water matric potential (m)
                       ! output: error control 
                       err,message)                                 ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 implicit none 
 ! input
 real(dp),intent(in)             :: stateVec(:)                     ! model state vector (mixed units)
 type(var_dlength),intent(in)    :: diag_data                       ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)    :: prog_data                       ! prognostic variables for a local HRU
 type(var_ilength),intent(in)    :: indx_data                       ! indices defining model states and layers                 
 ! output: variables for the vegetation canopy
 real(dp),intent(out)            :: scalarCanairTempTrial           ! trial value of canopy air temperature (K)
 real(dp),intent(out)            :: scalarCanopyTempTrial           ! trial value of canopy temperature (K)
 real(dp),intent(out)            :: scalarCanopyWatTrial            ! trial value of canopy total water (kg m-2)
 real(dp),intent(out)            :: scalarCanopyLiqTrial            ! trial value of canopy liquid water (kg m-2)
 real(dp),intent(out)            :: scalarCanopyIceTrial            ! trial value of canopy ice content (kg m-2)
 ! output: variables for the snow-soil domain
 real(dp),intent(out)            :: mLayerTempTrial(:)              ! trial vector of layer temperature (K)
 real(dp),intent(out)            :: mLayerVolFracWatTrial(:)        ! trial vector of volumetric total water content (-)
 real(dp),intent(out)            :: mLayerVolFracLiqTrial(:)        ! trial vector of volumetric liquid water content (-)
 real(dp),intent(out)            :: mLayerVolFracIceTrial(:)        ! trial vector of volumetric ice water content (-)
 real(dp),intent(out)            :: mLayerMatricHeadTrial(:)        ! trial vector of total water matric potential (m)
 real(dp),intent(out)            :: mLayerMatricHeadLiqTrial(:)     ! trial vector of liquid water matric potential (m)
 ! output: error control 
 integer(i4b),intent(out)        :: err                             ! error code
 character(*),intent(out)        :: message                         ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: iLayer                          ! index of layer within the snow+soil domain
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! number of model layers, and layer type
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
 ! indices defining model states and layers
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy energy state variable
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for energy states in the snow+soil subdomain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for hydrology states in the snow+soil subdomain
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):  [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):  [i4b]    number of hydrology variables in the snow+soil domain
 ! indices defining type of model state variables 
 ixStateType_subset      => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in):  [i4b(:)] [state subset] type of desired model state variables
 ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in):  [i4b(:)] index of the type of hydrology states in snow+soil domain
 ! model states for the vegetation canopy
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in):  [dp]     temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in):  [dp]     temperature of the vegetation canopy (K)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in):  [dp]     mass of total water on the vegetation canopy (kg m-2)
 ! model state variable vectors for the snow-soil layers
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in):  [dp(:)]  temperature of each snow/soil layer (K)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in):  [dp(:)]  volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in):  [dp(:)]  total water matric potential (m)
 mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,& ! intent(in):  [dp(:)]  liquid water matric potential (m)
 ! model diagnostic variables from a previous solution
 scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(in):  [dp(:)]  mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(in):  [dp(:)]  mass of ice on the vegetation canopy (kg m-2)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in):  [dp(:)]  volumetric fraction of liquid water (-)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat           & ! intent(in):  [dp(:)]  volumetric fraction of ice (-)
 ) ! association with variables in the data structures

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='varExtract/'

 ! *** extract state variables for the vegetation canopy

 ! initialize to state variable from the last update
 scalarCanairTempTrial = scalarCanairTemp
 scalarCanopyTempTrial = scalarCanopyTemp
 scalarCanopyWatTrial  = scalarCanopyWat
 scalarCanopyLiqTrial  = scalarCanopyLiq
 scalarCanopyIceTrial  = scalarCanopyIce

 ! check if computing the vegetation flux
 if(ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing .or. ixVegHyd/=integerMissing)then
 
  ! extract temperature of the canopy air space
  if(ixCasNrg/=integerMissing) scalarCanairTempTrial = stateVec(ixCasNrg)
  
  ! extract canopy temperature
  if(ixVegNrg/=integerMissing) scalarCanopyTempTrial = stateVec(ixVegNrg) 
  
  ! extract intercepted water
  if(ixVegHyd/=integerMissing)then
   select case( ixStateType_subset(ixVegHyd) )
    case(iname_liqCanopy); scalarCanopyLiqTrial = stateVec(ixVegHyd)
    case(iname_watCanopy); scalarCanopyWatTrial = stateVec(ixVegHyd)
    case default; err=20; message=trim(message)//'case not found: expect iname_liqCanopy or iname_watCanopy'; return
   end select
  endif
 
 endif  ! not computing the vegetation flux

 ! *** extract state variables from the snow+soil sub-domain

 ! initialize to the state variable from the last update
 mLayerTempTrial          = mLayerTemp
 mLayerVolFracWatTrial    = mLayerVolFracWat
 mLayerVolFracLiqTrial    = mLayerVolFracLiq
 mLayerVolFracIceTrial    = mLayerVolFracIce
 mLayerMatricHeadTrial    = mLayerMatricHead      ! total water matric potential
 mLayerMatricHeadLiqTrial = mLayerMatricHeadLiq   ! liquid water matric potential

 ! overwrite with the energy values from the state vector
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   mLayerTempTrial(iLayer) = stateVec( ixSnowSoilNrg(iLayer) ) 
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! overwrite with the hydrology values from the state vector
 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
   select case( ixHydType(iLayer) )
    case(iname_watLayer); mLayerVolFracWatTrial(iLayer)          = stateVec( ixSnowSoilHyd(iLayer) ) ! total water state variable for snow+soil layers
    case(iname_liqLayer); mLayerVolFracLiqTrial(iLayer)          = stateVec( ixSnowSoilHyd(iLayer) ) ! liquid water state variable for snow+soil layers
    case(iname_matLayer); mLayerMatricHeadTrial(iLayer-nSnow)    = stateVec( ixSnowSoilHyd(iLayer) ) ! total water matric potential variable for soil layers
    case(iname_lmpLayer); mLayerMatricHeadLiqTrial(iLayer-nSnow) = stateVec( ixSnowSoilHyd(iLayer) ) ! liquid matric potential state variable for soil layers
    case default ! do nothing
   end select
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 end associate

 end subroutine varExtract

end module getVectorz_module
