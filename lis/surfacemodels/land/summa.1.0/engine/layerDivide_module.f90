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

module layerDivide_module

! variable types
USE nrtype

! physical constants
USE multiconst,only:&
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)

! access named variables for snow and soil
USE globalData,only:iname_snow        ! named variables for snow
USE globalData,only:iname_soil        ! named variables for soil

! define look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers,       & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex       ! CLM option: combination/sub-dividion rules depend on layer index

! define look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:&
 noah_mp,                  & ! full Noah-MP implementation (including albedo)
 CLM_2stream,              & ! CLM 2-stream model (see CLM documentation)
 UEB_2stream,              & ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
 NL_scatter,               & ! Simplified method Nijssen and Lettenmaier (JGR 1999)
 BeersLaw                    ! Beer's Law (as implemented in VIC)

! define look-up values for the choice of albedo method
USE mDecisions_module,only:& ! identify model options for snow albedo
 constantDecay,            & ! constant decay in snow albedo (e.g., VIC, CLASS)
 variableDecay               ! variable decay in snow albedo (e.g., BATS approach, with destructive metamorphism + soot content)

implicit none
private
public::layerDivide

! provide access to the number layers throughout the module
integer(i4b)                    :: nSnow               ! number of snow layers
integer(i4b)                    :: nSoil               ! number of soil layers
integer(i4b)                    :: nLayers             ! total number of layers
! define missing values
real(dp)              :: missingDouble=-9999._dp  ! missing value (double precision)
integer(i4b)          :: missingInteger=-9999     ! missing value (integer)

contains

 ! ***********************************************************************************************************
 ! public subroutine layerDivide: add new snowfall to the system, and increase number of snow layers if needed
 ! ***********************************************************************************************************
 subroutine layerDivide(&
                        ! input/output: model data structures
                        model_decisions,             & ! intent(in):    model decisions
                        mpar_data,                   & ! intent(in):    model parameters
                        indx_data,                   & ! intent(inout): type of each layer
                        prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                        diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                        flux_data,                   & ! intent(inout): model fluxes for a local HRU
                        ! output
                        divideLayer,                 & ! intent(out): flag to denote that a layer was divided
                        err,message)                   ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 ! access the derived types to define the data structures
 USE data_types,only:&
                     var_d,            & ! data vector (dp)
                     var_ilength,      & ! data vector with variable length dimension (i4b)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! access metadata
 USE globalData,only:prog_meta,diag_meta,flux_meta,indx_meta   ! metadata
 ! access named variables defining elements in the data structures
 USE var_lookup,only:iLookPROG,iLookDIAG,iLookFLUX,iLookINDEX  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                            ! named variables for elements of the decision structure
 USE var_lookup,only:iLookPARAM                                ! named variables for elements of the parameter structure
 ! computational modules
 USE snow_utils_module,only:fracliquid,templiquid              ! functions to compute temperature/liquid water
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 type(var_dlength),intent(in)    :: mpar_data           ! model parameters
 type(var_ilength),intent(inout) :: indx_data           ! type of each layer
 type(var_dlength),intent(inout) :: prog_data           ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data           ! model diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data           ! model flux variables
 ! output
 logical(lgt),intent(out)        :: divideLayer         ! flag to denote that a layer was divided
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! define local variables
 character(LEN=256)              :: cmessage            ! error message of downwind routine
 integer(i4b)                    :: iLayer              ! layer index
 integer(i4b)                    :: jLayer              ! layer index
 real(dp),dimension(4)           :: zmax_lower          ! lower value of maximum layer depth
 real(dp),dimension(4)           :: zmax_upper          ! upper value of maximum layer depth
 real(dp)                        :: zmaxCheck           ! value of zmax for a given snow layer
 integer(i4b)                    :: nCheck              ! number of layers to check to divide
 logical(lgt)                    :: createLayer         ! flag to indicate we are creating a new snow layer
 real(dp)                        :: depthOriginal       ! original layer depth before sub-division (m)
 real(dp),parameter              :: fracTop=0.5_dp      ! fraction of old layer used for the top layer
 real(dp)                        :: surfaceLayerSoilTemp  ! temperature of the top soil layer (K)
 real(dp)                        :: maxFrozenSnowTemp   ! maximum temperature when effectively all water is frozen (K)
 real(dp),parameter              :: unfrozenLiq=0.01_dp ! unfrozen liquid water used to compute maxFrozenSnowTemp (-)
 real(dp)                        :: volFracWater        ! volumetric fraction of total water, liquid and ice (-)
 real(dp)                        :: fracLiq             ! fraction of liquid water (-)
 integer(i4b),parameter          :: ixVisible=1         ! named variable to define index in array of visible part of the spectrum
 integer(i4b),parameter          :: ixNearIR=2          ! named variable to define index in array of near IR part of the spectrum
 real(dp),parameter              :: verySmall=1.e-10_dp ! a very small number (used for error checking)
 ! --------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="layerDivide/"
 ! --------------------------------------------------------------------------------------------------------
 ! associate variables in the data structures
 associate(&
 ! model decisions
 ix_snowLayers          => model_decisions(iLookDECISIONS%snowLayers)%iDecision, & ! decision for snow combination
 ! model parameters (compute layer temperature)
 fc_param               => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1),       & ! freezing curve parameter for snow (K-1)
 ! model parameters (new snow density)
 newSnowDenMin          => mpar_data%var(iLookPARAM%newSnowDenMin)%dat(1),       & ! minimum new snow density (kg m-3)
 newSnowDenMult         => mpar_data%var(iLookPARAM%newSnowDenMult)%dat(1),      & ! multiplier for new snow density (kg m-3)
 newSnowDenScal         => mpar_data%var(iLookPARAM%newSnowDenScal)%dat(1),      & ! scaling factor for new snow density (K)
 ! model parameters (control the depth of snow layers)
 zmax                   => mpar_data%var(iLookPARAM%zmax)%dat(1),                & ! maximum layer depth (m)
 zmaxLayer1_lower       => mpar_data%var(iLookPARAM%zmaxLayer1_lower)%dat(1),    & ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
 zmaxLayer2_lower       => mpar_data%var(iLookPARAM%zmaxLayer2_lower)%dat(1),    & ! maximum layer depth for the 2nd layer when only 2 layers (m)
 zmaxLayer3_lower       => mpar_data%var(iLookPARAM%zmaxLayer3_lower)%dat(1),    & ! maximum layer depth for the 3rd layer when only 3 layers (m)
 zmaxLayer4_lower       => mpar_data%var(iLookPARAM%zmaxLayer4_lower)%dat(1),    & ! maximum layer depth for the 4th layer when only 4 layers (m)
 zmaxLayer1_upper       => mpar_data%var(iLookPARAM%zmaxLayer1_upper)%dat(1),    & ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
 zmaxLayer2_upper       => mpar_data%var(iLookPARAM%zmaxLayer2_upper)%dat(1),    & ! maximum layer depth for the 2nd layer when > 2 layers (m)
 zmaxLayer3_upper       => mpar_data%var(iLookPARAM%zmaxLayer3_upper)%dat(1),    & ! maximum layer depth for the 3rd layer when > 3 layers (m)
 zmaxLayer4_upper       => mpar_data%var(iLookPARAM%zmaxLayer4_upper)%dat(1),    & ! maximum layer depth for the 4th layer when > 4 layers (m)
 ! diagnostic scalar variables
 scalarSnowfall         => flux_data%var(iLookFLUX%scalarSnowfall)%dat(1),       & ! snowfall flux (kg m-2 s-1)
 scalarSnowfallTemp     => diag_data%var(iLookDIAG%scalarSnowfallTemp)%dat(1),   & ! computed temperature of fresh snow (K)
 scalarSnowDepth        => prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),      & ! total snow depth (m)
 scalarSWE              => prog_data%var(iLookPROG%scalarSWE)%dat(1)             & ! SWE (kg m-2)
 )  ! end associate statement

 ! ---------------------------------------------------------------------------------------------------

 ! initialize flag to denote that a layer was divided
 divideLayer=.false.

 ! identify algorithmic control parameters to syb-divide and combine snow layers
 zmax_lower = (/zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower/)
 zmax_upper = (/zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper/)

 ! initialize the number of snow layers
 nSnow   = indx_data%var(iLookINDEX%nSnow)%dat(1)
 nSoil   = indx_data%var(iLookINDEX%nSoil)%dat(1)
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)

 ! ***** special case of no snow layers
 if(nSnow==0)then

  ! check if create the first snow layer
  select case(ix_snowLayers)
   case(sameRulesAllLayers);    createLayer = (scalarSnowDepth > zmax)
   case(rulesDependLayerIndex); createLayer = (scalarSnowDepth > zmaxLayer1_lower)
   case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
  end select ! (option to combine/sub-divide snow layers)

  ! ** create a new snow layer
  if(createLayer)then

   ! flag that the layers have changed
   divideLayer=.true.

   ! add a layer to all model variables
   iLayer=0 ! (layer to divide: 0 is the special case of "snow without a layer")
   call addModelLayer(prog_data,prog_meta,iLayer,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
   call addModelLayer(diag_data,diag_meta,iLayer,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
   call addModelLayer(flux_data,flux_meta,iLayer,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
   call addModelLayer(indx_data,indx_meta,iLayer,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

   ! associate local variables to the information in the data structures
   ! NOTE: need to do this here, since state vectors have just been modified
   associate(&
   ! coordinate variables
   mLayerDepth      => prog_data%var(iLookPROG%mLayerDepth)%dat        ,& ! depth of each layer (m)
   ! model state variables
   mLayerTemp       => prog_data%var(iLookPROG%mLayerTemp)%dat         ,& ! temperature of each layer (K)
   mLayerVolFracIce => prog_data%var(iLookPROG%mLayerVolFracIce)%dat   ,& ! volumetric fraction of ice in each layer (-)
   mLayerVolFracLiq => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat    & ! volumetric fraction of liquid water in each layer (-)
   ) ! (association of local variables to the information in the data structures)

   ! get the layer depth
   mLayerDepth(1) = scalarSnowDepth

   ! compute surface layer temperature
   surfaceLayerSoilTemp = mLayerTemp(2)    ! temperature of the top soil layer (K)
   maxFrozenSnowTemp    = templiquid(unfrozenLiq,fc_param)               ! snow temperature at fraction "unfrozenLiq" (K)
   mLayerTemp(1)        = min(maxFrozenSnowTemp,surfaceLayerSoilTemp)    ! snow temperature  (K)

   ! compute the fraction of liquid water associated with the layer temperature
   fracLiq      = fracliquid(mLayerTemp(1),fc_param)

   ! compute volumeteric fraction of liquid water and ice
   volFracWater = (scalarSWE/scalarSnowDepth)/iden_water  ! volumetric fraction of total water (liquid and ice)
   mLayerVolFracIce(1) = (1._dp - fracLiq)*volFracWater*(iden_water/iden_ice)   ! volumetric fraction of ice (-)
   mLayerVolFracLiq(1) =          fracLiq *volFracWater                         ! volumetric fraction of liquid water (-)

   ! end association with local variables to the information in the data structures)
   end associate

   ! initialize albedo
   ! NOTE: albedo is computed within the Noah-MP radiation routine
   if(model_decisions(iLookDECISIONS%canopySrad)%iDecision /= noah_mp)then
    select case(model_decisions(iLookDECISIONS%alb_method)%iDecision)
     ! (constant decay rate -- albedo the same for all spectral bands)
     case(constantDecay)
      prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(1)          = mpar_data%var(iLookPARAM%albedoMax)%dat(1)
      prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(:) = mpar_data%var(iLookPARAM%albedoMax)%dat(1)
     ! (variable decay rate)
     case(variableDecay)
      prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(ixVisible) = mpar_data%var(iLookPARAM%albedoMaxVisible)%dat(1)
      prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(ixNearIR)  = mpar_data%var(iLookPARAM%albedoMaxNearIR)%dat(1)
      prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(1)                  = (        mpar_data%var(iLookPARAM%Frad_vis)%dat(1))*mpar_data%var(iLookPARAM%albedoMaxVisible)%dat(1) + &
                                                                          (1._dp - mpar_data%var(iLookPARAM%Frad_vis)%dat(1))*mpar_data%var(iLookPARAM%albedoMaxNearIR)%dat(1)
     case default; err=20; message=trim(message)//'unable to identify option for snow albedo'; return
    end select  ! identify option for snow albedo
    ! set direct albedo to diffuse albedo
    diag_data%var(iLookDIAG%spectralSnowAlbedoDirect)%dat(:) = prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(:)
   end if  ! (if NOT using the Noah-MP radiation routine)

  end if  ! if creating a new layer

 ! end special case of nSnow=0
 ! ********************************************************************************************************************
 ! ********************************************************************************************************************

 ! ***** sub-divide snow layers, if necessary
 else ! if nSnow>0

  ! identify the number of layers to check for need for sub-division
  select case(ix_snowLayers)
   case(sameRulesAllLayers);    nCheck = nSnow
   case(rulesDependLayerIndex); nCheck = min(nSnow,4)  ! the depth of the 5th layer, if it exists, does not have a maximum value
   case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
  end select ! (option to combine/sub-divide snow layers)
  
  ! loop through all layers, and sub-divide a given layer, if necessary
  do iLayer=1,nCheck
  
   ! identify the maximum depth of the layer
   select case(ix_snowLayers)
    case(sameRulesAllLayers);    zmaxCheck = zmax
    case(rulesDependLayerIndex)
     if(iLayer == nSnow)then
      zmaxCheck = zmax_lower(iLayer)
     else
      zmaxCheck = zmax_upper(iLayer)
     end if
    case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
   end select ! (option to combine/sub-divide snow layers)
  
   ! check the need to sub-divide
   if(prog_data%var(iLookPROG%mLayerDepth)%dat(iLayer) > zmaxCheck)then
  
    ! flag that layers were divided
    divideLayer=.true.
  
    ! add a layer to all model variables
    call addModelLayer(prog_data,prog_meta,iLayer,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    call addModelLayer(diag_data,diag_meta,iLayer,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    call addModelLayer(flux_data,flux_meta,iLayer,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    call addModelLayer(indx_data,indx_meta,iLayer,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
    ! define the layer depth
    layerSplit: associate(mLayerDepth => prog_data%var(iLookPROG%mLayerDepth)%dat)
    depthOriginal = mLayerDepth(iLayer)
    mLayerDepth(iLayer)   = fracTop*depthOriginal
    mLayerDepth(iLayer+1) = (1._dp - fracTop)*depthOriginal
    end associate layerSplit
  
    exit  ! NOTE: only sub-divide one layer per substep
  
   end if   ! (if sub-dividing layer)
  
  end do  ! (looping through layers)

 end if  ! if nSnow==0

 ! update coordinates
 if(divideLayer)then

  ! associate coordinate variables in data structure
  geometry: associate(&
  mLayerDepth      => prog_data%var(iLookPROG%mLayerDepth)%dat        ,& ! depth of the layer (m)
  mLayerHeight     => prog_data%var(iLookPROG%mLayerHeight)%dat       ,& ! height of the layer mid-point (m)
  iLayerHeight     => prog_data%var(iLookPROG%iLayerHeight)%dat       ,& ! height of the layer interface (m)
  layerType        => indx_data%var(iLookINDEX%layerType)%dat         ,& ! type of each layer (iname_snow or iname_soil)
  nSnow            => indx_data%var(iLookINDEX%nSnow)%dat(1)          ,& ! number of snow layers
  nSoil            => indx_data%var(iLookINDEX%nSoil)%dat(1)          ,& ! number of soil layers
  nLayers          => indx_data%var(iLookINDEX%nLayers)%dat(1)         & ! total number of layers
  )  ! (association of local variables with coordinate variab;es in data structures)

  ! update the layer type
  layerType(1:nSnow+1)         = iname_snow
  layerType(nSnow+2:nLayers+1) = iname_soil

  ! identify the number of snow and soil layers, and check all is a-OK
  nSnow   = count(layerType==iname_snow)
  nSoil   = count(layerType==iname_soil)
  nLayers = nSnow + nSoil

  ! re-set coordinate variables
  iLayerHeight(0) = -scalarSnowDepth
  do jLayer=1,nLayers
   iLayerHeight(jLayer) = iLayerHeight(jLayer-1) + mLayerDepth(jLayer)
   mLayerHeight(jLayer) = (iLayerHeight(jLayer-1) + iLayerHeight(jLayer))/2._dp
  end do

  ! check
  if(abs(sum(mLayerDepth(1:nSnow)) - scalarSnowDepth) > verySmall)then
   print*, 'nSnow = ', nSnow
   write(*,'(a,1x,f30.25,1x)') 'sum(mLayerDepth(1:nSnow)) = ', sum(mLayerDepth(1:nSnow))
   write(*,'(a,1x,f30.25,1x)') 'scalarSnowDepth           = ', scalarSnowDepth
   write(*,'(a,1x,f30.25,1x)') 'epsilon(scalarSnowDepth)  = ', epsilon(scalarSnowDepth)
   message=trim(message)//'sum of layer depths does not equal snow depth'
   err=20; return
  end if

  ! end association with coordinate variables in data structure
  end associate geometry

 end if  ! if dividing a layer

 ! end associate variables in data structure
 end associate

 end subroutine layerDivide


 ! ************************************************************************************************
 ! private subroutine addModelLayer: add an additional layer to all model vectors
 ! ************************************************************************************************
 subroutine addModelLayer(dataStruct,metaStruct,ix_divide,err,message)
 USE var_lookup,only:iLookVarType                  ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName        ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc             ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength       ! data vectors with variable length dimension
 USE data_types,only:var_info                      ! metadata structure
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 class(*),intent(inout)          :: dataStruct     ! data structure
 type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in)         :: ix_divide      ! index of the layer to divide
 ! output: error control
 integer(i4b),intent(out)        :: err            ! error code
 character(*),intent(out)        :: message        ! error message
 ! ---------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: ivar           ! index of model variable
 integer(i4b)                    :: ix_lower       ! lower bound of the vector
 integer(i4b)                    :: ix_upper       ! upper bound of the vector
 logical(lgt)                    :: stateVariable  ! .true. if variable is a state variable
 real(dp),allocatable            :: tempVec_dp(:)  ! temporary vector (double precision)
 integer(i4b),allocatable        :: tempVec_i4b(:) ! temporary vector (integer)
 character(LEN=256)              :: cmessage       ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='addModelLayer/'

 ! ***** add a layer to each model variable
 do ivar=1,size(metaStruct)
  
  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=nSnow
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=nSnow
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nLayers
   case default; cycle
  end select
  
  ! identify whether it is a state variable
  select case(trim(metaStruct(ivar)%varname))
   case('mLayerDepth','mLayerTemp','mLayerVolFracIce','mLayerVolFracLiq'); stateVariable=.true.
   case default; stateVariable=.false.
  end select

  ! divide layers
  select type(dataStruct)

   ! ** double precision
   type is (var_dlength)
    ! check allocated
    if(.not.allocated(dataStruct%var(ivar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! assign the data vector to the temporary vector
    call cloneStruc(tempVec_dp, ix_lower, source=dataStruct%var(ivar)%dat, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! reallocate space for the new vector
    deallocate(dataStruct%var(ivar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; end if
    allocate(dataStruct%var(ivar)%dat(ix_lower:ix_upper+1),stat=err)
    if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; end if
    ! populate the state vector
    if(stateVariable)then
     if(ix_upper > 0)then  ! (only copy data if the vector exists -- can be a variable for snow, with no layers)
      if(ix_divide > 0)then
       dataStruct%var(ivar)%dat(1:ix_divide)            = tempVec_dp(1:ix_divide)  ! copy data
       dataStruct%var(ivar)%dat(ix_divide+1)            = tempVec_dp(ix_divide)    ! repeat data for the sub-divided layer
      end if
      if(ix_upper > ix_divide) &
       dataStruct%var(ivar)%dat(ix_divide+2:ix_upper+1) = tempVec_dp(ix_divide+1:ix_upper)  ! copy data
     end if  ! if the vector exists
    ! not a state variable
    else
     dataStruct%var(ivar)%dat(:) = missingDouble
    end if
    ! deallocate the temporary vector: strictly not necessary, but include to be safe
    deallocate(tempVec_dp,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! ** integer
   type is (var_ilength)
    ! check allocated
    if(.not.allocated(dataStruct%var(ivar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! assign the data vector to the temporary vector
    call cloneStruc(tempVec_i4b, ix_lower, source=dataStruct%var(ivar)%dat, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! reallocate space for the new vector
    deallocate(dataStruct%var(ivar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; end if
    allocate(dataStruct%var(ivar)%dat(ix_lower:ix_upper+1),stat=err)
    if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; end if
    ! populate the state vector
    if(stateVariable)then
     if(ix_upper > 0)then  ! (only copy data if the vector exists -- can be a variable for snow, with no layers)
      if(ix_divide > 0)then
       dataStruct%var(ivar)%dat(1:ix_divide)            = tempVec_i4b(1:ix_divide)  ! copy data
       dataStruct%var(ivar)%dat(ix_divide+1)            = tempVec_i4b(ix_divide)    ! repeat data for the sub-divided layer
      end if
      if(ix_upper > ix_divide) &
       dataStruct%var(ivar)%dat(ix_divide+2:ix_upper+1) = tempVec_i4b(ix_divide+1:ix_upper)  ! copy data
     end if  ! if the vector exists
    ! not a state variable
    else
     dataStruct%var(ivar)%dat(:) = missingInteger
    end if
    ! deallocate the temporary vector: strictly not necessary, but include to be safe
    deallocate(tempVec_i4b,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! check that we found the data type
   class default; err=20; message=trim(message)//'unable to identify the data type'; return

  end select ! dependence on data types

 end do  ! looping through variables

 end subroutine addModelLayer

end module layerDivide_module
