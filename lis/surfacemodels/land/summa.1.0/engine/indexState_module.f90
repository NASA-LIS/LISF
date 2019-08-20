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

module indexState_module

! data types
USE nrtype

! missing data
USE globalData,only:integerMissing  ! missing integer

! named variables for domain types
USE globalData,only:iname_cas       ! canopy air space
USE globalData,only:iname_veg       ! vegetation
USE globalData,only:iname_snow      ! snow
USE globalData,only:iname_soil      ! soil

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

! provide access to the derived types to define the data structures
USE data_types,only:var_ilength     ! data vector with variable length dimension (i4b)

! provide access to the metadata
USE globalData,only:indx_meta       ! metadata for the variables in the index structure

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! provide access to the missing f2008 functions
USE f2008funcs_module,only:findIndex              ! finds the index of the first value within a vector

! provide access to the numerical recipes utility modules
USE nr_utility_module,only:arth                           ! creates a sequence of numbers (start, incr, n)

implicit none
private
public::indexState
public::indexSplit
contains


 ! **********************************************************************************************************
 ! public subroutine indexState: define list of indices for each state variable 
 ! **********************************************************************************************************
 subroutine indexState(computeVegFlux,          & ! intent(in):    flag to denote if computing the vegetation flux
                       nSnow,nSoil,nLayers,     & ! intent(in):    number of snow and soil layers, and total number of layers
                       indx_data,               & ! intent(inout): indices defining model states and layers
                       err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input
 logical(lgt),intent(in)         :: computeVegFlux         ! flag to denote if computing the vegetation flux
 integer(i4b),intent(in)         :: nSnow,nSoil,nLayers    ! number of snow and soil layers, and total number of layers
 type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! general local variables
 character(len=256)              :: cmessage               ! message of downwind routine
 integer(i4b),parameter          :: nVarSnowSoil=2         ! number of state variables in the snow and soil domain (energy and total water/matric head)
 ! indices of model state variables
 integer(i4b)                    :: ixTopNrg               ! index of upper-most energy state in the snow-soil subdomain
 integer(i4b)                    :: ixTopWat               ! index of upper-most total water state in the snow-soil subdomain
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! number of state variables of different type
 nCasNrg       => indx_data%var(iLookINDEX%nVegNrg)%dat(1)   , & ! number of energy state variables for the canopy air space
 nVegNrg       => indx_data%var(iLookINDEX%nVegNrg)%dat(1)   , & ! number of energy state variables for the vegetation canopy
 nVegMass      => indx_data%var(iLookINDEX%nVegMass)%dat(1)  , & ! number of hydrology states for vegetation (mass of water)
 nVegState     => indx_data%var(iLookINDEX%nVegState)%dat(1) , & ! number of vegetation state variables
 nNrgState     => indx_data%var(iLookINDEX%nNrgState)%dat(1) , & ! number of energy state variables
 nWatState     => indx_data%var(iLookINDEX%nWatState)%dat(1) , & ! number of "total water" states (vol. total water content)
 nMatState     => indx_data%var(iLookINDEX%nMatState)%dat(1) , & ! number of matric head state variables
 nMassState    => indx_data%var(iLookINDEX%nMassState)%dat(1), & ! number of hydrology state variables (mass of water)
 nState        => indx_data%var(iLookINDEX%nState)%dat(1)    , & ! total number of model state variables
 ! vectors of indices for specfic state types within specific sub-domains IN THE FULL STATE VECTOR
 ixNrgCanair   => indx_data%var(iLookINDEX%ixNrgCanair)%dat  , & ! indices IN THE FULL VECTOR for energy states in canopy air space domain
 ixNrgCanopy   => indx_data%var(iLookINDEX%ixNrgCanopy)%dat  , & ! indices IN THE FULL VECTOR for energy states in the canopy domain
 ixHydCanopy   => indx_data%var(iLookINDEX%ixHydCanopy)%dat  , & ! indices IN THE FULL VECTOR for hydrology states in the canopy domain
 ixNrgLayer    => indx_data%var(iLookINDEX%ixNrgLayer)%dat   , & ! indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer    => indx_data%var(iLookINDEX%ixHydLayer)%dat   , & ! indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ! indices for model state variables
 ixSoilState   => indx_data%var(iLookINDEX%ixSoilState)%dat  , & ! list of indices for all soil layers
 ixLayerState  => indx_data%var(iLookINDEX%ixLayerState)%dat   & ! list of indices for all model layers
 ) ! association to variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='indexState/'

 ! -----
 ! * define the number of state variables...
 ! -----------------------------------------

 ! define the number of vegetation state variables (defines position of snow-soil states in the state vector)
 if(computeVegFlux)then
  nCasNrg   = 1
  nVegNrg   = 1
  nVegMass  = 1
  nVegState = nCasNrg + nVegNrg + nVegMass
 else
  nCasNrg   = 0
  nVegNrg   = 0
  nVegMass  = 0
  nVegState = 0
 end if

 ! define the number state variables of different type
 nNrgState  = nCasNrg + nVegNrg + nLayers  ! number of energy state variables
 nWatState  = nSnow                        ! number of "total water" state variables -- will be modified later if using primary variable switching 
 nMatState  = nSoil                        ! number of matric head state variables -- will be modified later if using primary variable switching
 nMassState = nVegMass                     ! number of mass state variables -- currently restricted to canopy water

 ! define the number of model state variables
 nState = nVegState + nLayers*nVarSnowSoil   ! *nVarSnowSoil (both energy and total water)

 ! -----
 ! * define the indices of state variables WITHIN THE FULL STATE VECTOR...
 ! -----------------------------------------------------------------------

 ! define indices in the vegetation domain
 if(computeVegFlux)then
  ixNrgCanair = 1 ! indices IN THE FULL VECTOR for energy states in canopy air space domain  (-)
  ixNrgCanopy = 2 ! indices IN THE FULL VECTOR for energy states in the canopy domain        (-)
  ixHydCanopy = 3 ! indices IN THE FULL VECTOR for hydrology states in the canopy domain     (-)
 else
  ixNrgCanair = integerMissing
  ixNrgCanopy = integerMissing
  ixHydCanopy = integerMissing
 end if

 ! define the index of the top layer
 ! NOTE: local variables -- actual indices defined when building the state subset
 ixTopNrg = nVegState + 1                       ! energy
 ixTopWat = nVegState + 2                       ! total water (only snow)

 ! define the indices within the snow+soil domain
 ixNrgLayer = arth(ixTopNrg,nVarSnowSoil,nLayers)  ! energy
 ixHydLayer = arth(ixTopWat,nVarSnowSoil,nLayers)  ! total water

 ! -----
 ! * define the type of model states...
 ! ------------------------------------

 ! re-allocate index vectors for the full state vector (if needed)...
 call resizeIndx( (/iLookINDEX%ixMapFull2Subset, iLookINDEX%ixControlVolume, iLookINDEX%ixDomainType, iLookINDEX%ixStateType, iLookINDEX%ixAllState/), & ! desired variables
                  indx_data,  & ! data structure
                  nState,     & ! vector length
                  err,cmessage) ! error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! make an association to the ALLOCATABLE variables in the data structures
 ! NOTE: we need to do this here since the size may have changed above
 associate(&
 ixControlVolume => indx_data%var(iLookINDEX%ixControlVolume)%dat , & ! index of control volume for different domains (veg, snow, soil)
 ixDomainType    => indx_data%var(iLookINDEX%ixDomainType)%dat    , & ! indices defining the type of the domain (iname_veg, iname_snow, iname_soil)
 ixStateType     => indx_data%var(iLookINDEX%ixStateType)%dat     , & ! indices defining the type of the state (iname_nrgLayer...)
 ixAllState      => indx_data%var(iLookINDEX%ixAllState)%dat        & ! list of indices for all model state variables
 )  ! making an association to variables in the data structures

 ! define indices for state variables
 ixAllState   = arth(1,1,nState)
 ixSoilState  = arth(1,1,nSoil)
 ixLayerState = arth(1,1,nLayers)

 ! define the state type for the vegetation canopy
 if(computeVegFlux)then
  ixStateType(ixNrgCanair) = iname_nrgCanair
  ixStateType(ixNrgCanopy) = iname_nrgCanopy
  ixStateType(ixHydCanopy) = iname_watCanopy
 endif

 ! define the state type for the snow+soil domain (energy)
 ixStateType(ixNrgLayer) = iname_nrgLayer

 ! define the state type for the snow+soil domain (hydrology)
 if(nSnow>0) ixStateType( ixHydLayer(      1:nSnow)   ) = iname_watLayer
             ixStateType( ixHydLayer(nSnow+1:nLayers) ) = iname_matLayer ! refine later to be either iname_watLayer or iname_matLayer

 ! define the domain type for vegetation
 if(computeVegFlux)then
  ixDomainType(ixNrgCanair) = iname_cas
  ixDomainType(ixNrgCanopy) = iname_veg
  ixDomainType(ixHydCanopy) = iname_veg
 endif

 ! define the domain type for snow
 if(nSnow>0)then
  ixDomainType( ixNrgLayer(1:nSnow) ) = iname_snow
  ixDomainType( ixHydLayer(1:nSnow) ) = iname_snow
 endif

 ! define the domain type for soil
 ixDomainType( ixNrgLayer(nSnow+1:nLayers) ) = iname_soil
 ixDomainType( ixHydLayer(nSnow+1:nLayers) ) = iname_soil

 ! define the index of each control volume in the vegetation domains
 if(computeVegFlux)then
  ixControlVolume(ixNrgCanair) = 1  ! NOTE: assumes scalar
  ixControlVolume(ixNrgCanopy) = 1
  ixControlVolume(ixHydCanopy) = 1
 endif

 ! define the index of the each control volume in the snow domain
 if(nSnow>0)then
  ixControlVolume( ixNrgLayer(1:nSnow) ) = ixLayerState(1:nSnow)
  ixControlVolume( ixHydLayer(1:nSnow) ) = ixLayerState(1:nSnow)
 endif

 ! define the index of the each control volume in the soil domain
 ixControlVolume( ixNrgLayer(nSnow+1:nLayers) ) = ixSoilState(1:nSoil)
 ixControlVolume( ixHydLayer(nSnow+1:nLayers) ) = ixSoilState(1:nSoil)

 !print*, 'ixControlVolume = ', ixControlVolume
 !print*, 'ixDomainType    = ', ixDomainType
 !print*, 'ixStateType     = ', ixStateType

 ! end association to the ALLOCATABLE variables in the data structures
 end associate 

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 end associate  ! end association to variables in the data structures
 end subroutine indexState


 ! **********************************************************************************************************
 ! public subroutine indexSplit: define list of indices for each state variable 
 ! **********************************************************************************************************
 subroutine indexSplit(stateSubsetMask,             & ! intent(in)    : logical vector (.true. if state is in the subset)
                       nSnow,nSoil,nLayers,nSubset, & ! intent(in)    : number of snow and soil layers, and total number of layers
                       indx_data,                   & ! intent(inout) : index data structure
                       err,message)                   ! intent(out)   : error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input
 logical(lgt),intent(in)         :: stateSubsetMask(:)          ! logical vector (.true. if state is in the subset) 
 integer(i4b),intent(in)         :: nSnow,nSoil,nLayers,nSubset ! number of snow and soil layers, total number of layers, and number of states in the subset
 type(var_ilength),intent(inout) :: indx_data                   ! indices defining model states and layers
 ! output
 integer(i4b),intent(out)        :: err                         ! error code
 character(*),intent(out)        :: message                     ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: iVar                        ! variable index
 integer(i4b)                    :: ixVegWat                    ! index of total water in the vegetation canopy
 integer(i4b)                    :: ixVegLiq                    ! index of liquid water in the vegetation canopy
 integer(i4b)                    :: ixTopWat                    ! index of upper-most total water state in the snow-soil subdomain
 integer(i4b)                    :: ixTopLiq                    ! index of upper-most liquid water state in the snow-soil subdomain
 integer(i4b)                    :: ixTopMat                    ! index of upper-most total water matric potential state in the soil subdomain
 integer(i4b)                    :: ixTopLMP                    ! index of upper-most liquid water matric potential state in the soil subdomain
 integer(i4b),dimension(nSubset) :: ixSequence                  ! sequential index in model state vector
 logical(lgt),dimension(nSubset) :: stateTypeMask               ! mask of state vector for specific state subsets
 logical(lgt),dimension(nLayers) :: volFracWat_mask             ! mask of layers within the snow+soil domain
 logical(lgt),dimension(nSoil)   :: matricHead_mask             ! mask of layers within the soil domain
 character(len=256)              :: cmessage                    ! error message of downwind routine
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! make association to variables in the data structures
 fullState: associate(&

 ! indices of model state variables for the vegetation domain
 ixCasNrg         => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)      ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
 ixVegNrg         => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)      ,& ! intent(in):    [i4b]    index of canopy energy state variable
 ixVegHyd         => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)      ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)

 ! indices of the top model state variables in the snow+soil system
 ixTopNrg         => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)      ,& ! intent(in):    [i4b]    index of upper-most energy state in the snow-soil subdomain
 ixTopHyd         => indx_data%var(iLookINDEX%ixTopHyd)%dat(1)      ,& ! intent(in):    [i4b]    index of upper-most hydrology state in the snow-soil subdomain
 
 ! indices of model state variables
 ixMapFull2Subset => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat ,& ! intent(in):    [i4b(:)] list of indices in the state subset (missing for values not in the subset)
 ixDomainType     => indx_data%var(iLookINDEX%ixDomainType)%dat     ,& ! intent(in):    [i4b(:)] indices defining the domain of the state (iname_veg, iname_snow, iname_soil)
 ixStateType      => indx_data%var(iLookINDEX%ixStateType)%dat      ,& ! intent(in):    [i4b(:)] indices defining the type of the state (ixNrgState...)
 ixAllState       => indx_data%var(iLookINDEX%ixAllState)%dat       ,& ! intent(in):    [i4b(:)] list of indices for all model state variables (1,2,3,...nState)
 ixNrgLayer       => indx_data%var(iLookINDEX%ixNrgLayer)%dat       ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer       => indx_data%var(iLookINDEX%ixHydLayer)%dat       ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ixHydType        => indx_data%var(iLookINDEX%ixHydType)%dat        ,& ! intent(in):    [i4b(:)] index of the type of hydrology states in snow+soil domain

 ! indices of the entire state vector, all model layers, and soil layers
 ixSoilState      => indx_data%var(iLookINDEX%ixSoilState)%dat      ,& ! intent(in):    [i4b(:)] list of indices for all soil layers
 ixLayerState     => indx_data%var(iLookINDEX%ixLayerState)%dat     ,& ! intent(in):    [i4b(:)] list of indices for all model layers

 ! vector of energy indices for the snow and soil domains
 ! NOTE: states not in the subset are equal to integerMissing
 ixSnowSoilNrg    => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat    ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowOnlyNrg    => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat    ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow domain
 ixSoilOnlyNrg    => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat    ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the soil domain

 ! vector of hydrology indices for the snow and soil domains
 ! NOTE: states not in the subset are equal to integerMissing
 ixSnowSoilHyd    => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat    ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 ixSnowOnlyHyd    => indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat    ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow domain
 ixSoilOnlyHyd    => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat    ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the soil domain

 ! number of state variables of a specific type
 nSnowSoilNrg     => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1) ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
 nSnowOnlyNrg     => indx_data%var(iLookINDEX%nSnowOnlyNrg )%dat(1) ,& ! intent(in):    [i4b]    number of energy state variables in the snow domain
 nSoilOnlyNrg     => indx_data%var(iLookINDEX%nSoilOnlyNrg )%dat(1) ,& ! intent(in):    [i4b]    number of energy state variables in the soil domain
 nSnowSoilHyd     => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1) ,& ! intent(in):    [i4b]    number of hydrology variables in the snow+soil domain
 nSnowOnlyHyd     => indx_data%var(iLookINDEX%nSnowOnlyHyd )%dat(1) ,& ! intent(in):    [i4b]    number of hydrology variables in the snow domain
 nSoilOnlyHyd     => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)  & ! intent(in):    [i4b]    number of hydrology variables in the soil domain
 
 ) ! association to variables in the data structures

 ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='indexSplit/'

 ! -----
 ! - preliminaries...
 ! ------------------

 ! define the type of variable in the snow+soil domain
 ixHydType(1:nLayers) = ixStateType( ixHydLayer(1:nLayers) )

 ! get the mapping between the full state vector and the state subset
 ixMapFull2Subset( pack(ixAllState,      stateSubsetMask) ) = arth(1,1,nSubset)  ! indices in the state subset
 ixMapFull2Subset( pack(ixAllState, .not.stateSubsetMask) ) = integerMissing

 ! -----
 ! - get vectors of different state subsets...
 ! -------------------------------------------

 ! get different masks
 volFracWat_mask = (ixHydType==iname_watLayer .or. ixHydType==iname_liqLayer)
 matricHead_mask = (ixHydType(nSnow+1:nLayers)==iname_matLayer .or. ixHydType(nSnow+1:nLayers)==iname_lmpLayer)

 ! get state subsets for desired variables
 do iVar=1,size(indx_data%var)   ! loop through index variables

  ! get the subset of indices
  ! NOTE: indxSubset(subset, fullVector, mask), provides subset of fullVector where mask==.true.
  select case(iVar)
   case(iLookINDEX%ixMapSubset2Full);     call indxSubset(indx_data%var(iVar)%dat, ixAllState,   stateSubsetMask, err, cmessage) 
   case(iLookINDEX%ixStateType_subset);   call indxSubset(indx_data%var(iVar)%dat, ixStateType,  stateSubsetMask, err, cmessage)
   case(iLookINDEX%ixDomainType_subset);  call indxSubset(indx_data%var(iVar)%dat, ixDomainType, stateSubsetMask, err, cmessage)
   case(iLookINDEX%ixVolFracWat);         call indxSubset(indx_data%var(iVar)%dat, ixLayerState, volFracWat_mask, err, cmessage)
   case(iLookINDEX%ixMatricHead);         call indxSubset(indx_data%var(iVar)%dat, ixSoilState,  matricHead_mask, err, cmessage)
   case default; cycle ! only need to process the above variables
  end select  ! iVar
  if(err/=0)then; message=trim(message)//trim(cmessage)//'[varname='//trim(indx_meta(ivar)%varname)//']'; return; endif

 end do  ! looping through variables in the data structure

 ! make association to variables in the data structures
 subsetState: associate(ixStateType_subset => indx_data%var(iLookINDEX%ixStateType_subset)%dat) ! named variables defining the states in the subset

 ! -----
 ! - get indices for the (currently) scalar states in the vegetation domain...
 ! ---------------------------------------------------------------------------

 ! check the number of state variables in the vegetation canopy
 if(count(ixStateType_subset==iname_nrgCanair)>1)then; err=20; message=trim(message)//'expect count(iname_nrgCanair)=1 or 0'; return; endif  
 if(count(ixStateType_subset==iname_nrgCanopy)>1)then; err=20; message=trim(message)//'expect count(iname_nrgCanopy)=1 or 0'; return; endif  
 if(count(ixStateType_subset==iname_watCanopy)>1)then; err=20; message=trim(message)//'expect count(iname_watCanopy)=1 or 0'; return; endif  

 ! define indices for energy states for the canopy air space and the vegetation canopy
 ! NOTE: finds first index of named variable within stateType (set to integerMissing if not found)
 ixCasNrg = findIndex(ixStateType_subset, iname_nrgCanair, integerMissing)   ! energy of the canopy air space
 ixVegNrg = findIndex(ixStateType_subset, iname_nrgCanopy, integerMissing)   ! energy of the vegetation canopy

 ! define indices for hydrology states for the vegetation canopy
 ! NOTE: local variables -- ixVegHyd defined next
 ixVegWat = findIndex(ixStateType_subset, iname_watCanopy, integerMissing)   ! total water in the vegetation canopy
 ixVegLiq = findIndex(ixStateType_subset, iname_liqCanopy, integerMissing)   ! liquid water in the vegetation canopy
 ixVegHyd = merge(ixVegWat, ixVegLiq, ixVegWat/=integerMissing)

 ! define index for the upper-most energy state variables in the snow+soil domain
 ixTopNrg = findIndex(ixStateType_subset, iname_nrgLayer, integerMissing)    ! upper-most energy state in the snow+soil system

 ! define index for the upper-most hydrology state variables in the snow+soil domain
 ! NOTE: local variables -- ixTopHyd defined next
 ixTopWat = findIndex(ixStateType_subset, iname_watLayer, integerMissing)    ! upper-most total water state variable in the snow+soil system
 ixTopLiq = findIndex(ixStateType_subset, iname_liqLayer, integerMissing)    ! upper-most liquid water state variable in the snow+soil system
 ixTopMat = findIndex(ixStateType_subset, iname_matLayer, integerMissing)    ! upper-most total water matric potential state 
 ixTopLMP = findIndex(ixStateType_subset, iname_lmpLayer, integerMissing)    ! upper-most liquid water matric potential state 

 ! define index for the upper most hydrology state in the snow+soil system
 if(ixTopWat==integerMissing .and. ixTopLiq==integerMissing)then
  ixTopHyd = merge(ixTopMat, ixTopLMP, ixTopMat/=integerMissing)      ! no water state, so upper-most hydrology state is the upper-most matric head state (if it exists)
 else
  ixTopHyd = merge(ixTopWat, ixTopLiq, ixTopWat/=integerMissing)      ! ixTopWat is used if it is not missing
 endif

 ! -----
 ! - get vector of indices within the state subset state variables of a given type...
 ! ----------------------------------------------------------------------------------

 ! define index in full state vector
 ixSequence = arth(1,1,nSubset)

 ! get state subsets for desired variables
 do iVar=1,size(indx_data%var)   ! loop through index variables

  ! define the mask
  select case(iVar)
   case(iLookINDEX%ixNrgOnly);     stateTypeMask = (ixStateType_subset==iname_nrgCanair .or. ixStateType_subset==iname_nrgCanopy .or. ixStateType_subset==iname_nrgLayer)  ! list of indices for all energy states
   case(iLookINDEX%ixHydOnly);     stateTypeMask = (ixStateType_subset==iname_watLayer  .or. ixStateType_subset==iname_liqLayer  .or. ixStateType_subset==iname_matLayer .or. ixStateType_subset==iname_lmpLayer)   ! list of indices for all hydrology states
   case(iLookINDEX%ixMatOnly);     stateTypeMask = (ixStateType_subset==iname_matLayer  .or. ixStateType_subset==iname_lmpLayer)   ! list of indices for matric head state variables
   case(iLookINDEX%ixMassOnly);    stateTypeMask = (ixStateType_subset==iname_watCanopy)   ! list of indices for hydrology states (mass of water)
   case default; cycle ! only need to process the above variables
  end select  ! iVar

  ! get the subset of indices
  ! NOTE: indxSubset(subset, fullVector, mask), provides subset of fullVector where mask==.true.
  call indxSubset(indx_data%var(iVar)%dat,ixSequence,stateTypeMask,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage)//'[varname='//trim(indx_meta(ivar)%varname)//']'; return; endif

 end do  ! looping through variables in the data structure

 ! -----
 ! - get vector of indices of the state subset for layers in the snow+soil domain...
 ! ---------------------------------------------------------------------------------

 ! get list of indices for energy
 ! NOTE: layers not in the state subset will be missing
 ixSnowSoilNrg = ixMapFull2Subset(ixNrgLayer)                    ! both snow and soil layers
 ixSnowOnlyNrg = ixMapFull2Subset(ixNrgLayer(      1:nSnow  ))   ! snow layers only
 ixSoilOnlyNrg = ixMapFull2Subset(ixNrgLayer(nSnow+1:nLayers))   ! soil layers only

 ! get list of indices for hydrology
 ! NOTE: layers not in the state subset will be missing
 ixSnowSoilHyd = ixMapFull2Subset(ixHydLayer)                    ! both snow and soil layers
 ixSnowOnlyHyd = ixMapFull2Subset(ixHydLayer(      1:nSnow  ))   ! snow layers only
 ixSoilOnlyHyd = ixMapFull2Subset(ixHydLayer(nSnow+1:nLayers))   ! soil layers only

 ! get the number of valid states for energy
 nSnowSoilNrg = count(ixSnowSoilNrg/=integerMissing)
 nSnowOnlyNrg = count(ixSnowOnlyNrg/=integerMissing)
 nSoilOnlyNrg = count(ixSoilOnlyNrg/=integerMissing)

 ! get the number of valid states for hydrology
 nSnowSoilHyd = count(ixSnowSoilHyd/=integerMissing)
 nSnowOnlyHyd = count(ixSnowOnlyHyd/=integerMissing)
 nSoilOnlyHyd = count(ixSoilOnlyHyd/=integerMissing)

 ! end association to data in structures
 end associate subsetState
 end associate fullState

 end subroutine indexSplit


 ! **********************************************************************************************************
 ! private subroutine indxSubset: get a subset of indices for a given mask
 ! **********************************************************************************************************
 subroutine indxSubset(ixSubset,ixMaster,mask,err,message)
 implicit none
 ! input-output: subset of indices for allocation/population
 integer(i4b),intent(inout),allocatable :: ixSubset(:)           ! subset of indices
 ! input
 integer(i4b),intent(in)                :: ixMaster(:)           ! full list of indices
 logical(lgt),intent(in)                :: mask(:)               ! desired indices
 ! error control
 integer(i4b),intent(out)               :: err                   ! error code
 character(*),intent(out)               :: message               ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                           :: nSubset               ! length of the subset
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="indxSubset/"

 ! check size match
 if(size(ixMaster)/=size(mask))then
  message=trim(message)//'size mismatch'
  err=20; return
 endif

 ! get the number of variables
 nSubset = count(mask)

 ! check if we need to reallocate space
 if(size(ixSubset)/=nSubset) then

  ! deallocate space
  deallocate(ixSubset,stat=err)
  if(err/=0)then; message=trim(message)//'unable to deallocate space for variable'; err=20; return; endif

  ! allocate space
  allocate(ixSubset(nSubset),stat=err)
  if(err/=0)then; message=trim(message)//'unable to deallocate space for variable'; err=20; return; endif

 endif  ! allocating space

 ! define indices for variable types in specific sub-domains
 if(nSubset>0) ixSubset = pack(ixMaster, mask)

 end subroutine indxSubset





 ! **********************************************************************************************************
 ! private subroutine resizeIndx: re-size specific index vectors 
 ! **********************************************************************************************************
 subroutine resizeIndx(ixDesire,indx_data,nVec,err,message)
 ! input
 integer(i4b)     ,intent(in)    :: ixDesire(:)            ! variables needing to be re-sized
 type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
 integer(i4b)     ,intent(in)    :: nVec                   ! desired vector length 
 ! output
 integer(i4b)     ,intent(out)   :: err                    ! error code
 character(*)     ,intent(out)   :: message                ! error message 
 ! local variables
 integer(i4b)                    :: jVar,iVar              ! vatiable index
 ! initialize error control
 err=0; message='resizeIndx/'

 ! loop through variables
 do jVar=1,size(ixDesire)

  ! define index in index array
  iVar = ixDesire(jVar)

  ! check iVar is within range
  if(iVar<1 .or. iVar>size(indx_data%var))then
   message=trim(message)//'desired variable is out of range'
   err=20; return
  endif

  ! check if we need to reallocate space
  if(size(indx_data%var(iVar)%dat) == nVec) cycle

  ! deallocate space
  deallocate(indx_data%var(iVar)%dat,stat=err)
  if(err/=0)then
   message=trim(message)//'unable to deallocate space for variable '//trim(indx_meta(ivar)%varname)
   err=20; return
  endif

  ! allocate space
  allocate(indx_data%var(iVar)%dat(nVec),stat=err)
  if(err/=0)then
   message=trim(message)//'unable to allocate space for variable '//trim(indx_meta(ivar)%varname)
   err=20; return
  endif

  ! set to missing
  indx_data%var(iVar)%dat = integerMissing

 end do  ! looping through variables

 end subroutine resizeIndx

end module indexState_module
