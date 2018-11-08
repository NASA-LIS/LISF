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

module layerMerge_module
! data types
USE nrtype
! access named variables for snow and soil
USE globalData,only:iname_snow        ! named variables for snow
USE globalData,only:iname_soil        ! named variables for soil
! physical constants
USE multiconst,only:&
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)
! look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
! provide access to external modules
USE var_derive_module,only:calcHeight ! module to calculate height at layer interfaces and layer mid-point
implicit none
private
public::layerMerge
! provide access to the number layers throughout the module
integer(i4b)          :: nSnow                    ! number of snow layers
integer(i4b)          :: nSoil                    ! number of soil layers
integer(i4b)          :: nLayers                  ! total number of layers
! define missing values
real(dp)              :: missingDouble=-9999._dp  ! missing value (double precision)
integer(i4b)          :: missingInteger=-9999     ! missing value (integer)

contains


 ! *****************************************************************************************************************
 ! public subroutine layerMerge: merge layers if the thickness is less than zmin
 ! *****************************************************************************************************************
 subroutine layerMerge(&
                       ! input/output: model data structures
                       tooMuchMelt,                 & ! intent(in):    flag to force merge of snow layers
                       model_decisions,             & ! intent(in):    model decisions
                       mpar_data,                   & ! intent(in):    model parameters
                       indx_data,                   & ! intent(inout): type of each layer
                       prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,                   & ! intent(inout): model fluxes for a local HRU
                       ! output
                       mergedLayers,                & ! intent(out): flag to denote that layers were merged
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
 USE var_lookup,only:iLookPARAM,iLookPROG,iLookINDEX  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                   ! named variables for elements of the decision structure
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 logical(lgt),intent(in)         :: tooMuchMelt         ! flag to denote that ice is insufficient to support melt
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 type(var_dlength),intent(in)    :: mpar_data           ! model parameters
 type(var_ilength),intent(inout) :: indx_data           ! type of each layer
 type(var_dlength),intent(inout) :: prog_data           ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data           ! model diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data           ! model flux variables
 ! output
 logical(lgt),intent(out)        :: mergedLayers        ! flag to denote that layers were merged
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! define local variables
 character(LEN=256)              :: cmessage            ! error message of downwind routine
 real(dp),dimension(5)           :: zminLayer           ! minimum layer depth in each layer (m)
 logical(lgt)                    :: removeLayer         ! flag to indicate need to remove a layer
 integer(i4b)                    :: nCheck              ! number of layers to check for combination
 integer(i4b)                    :: iSnow               ! index of snow layers (looping)
 integer(i4b)                    :: jSnow               ! index of snow layer identified for combination with iSnow
 integer(i4b)                    :: kSnow               ! index of the upper layer of the two layers identified for combination
 ! initialize error control
 err=0; message="layerMerge/"
 ! --------------------------------------------------------------------------------------------------------
 ! associate variables to the data structures
 associate(&

 ! model decisions
 ix_snowLayers    => model_decisions(iLookDECISIONS%snowLayers)%iDecision, & ! decision for snow combination

 ! model parameters (control the depth of snow layers)
 zmin             => mpar_data%var(iLookPARAM%zmin)%dat(1),                & ! minimum layer depth (m)
 zminLayer1       => mpar_data%var(iLookPARAM%zminLayer1)%dat(1),          & ! minimum layer depth for the 1st (top) layer (m)
 zminLayer2       => mpar_data%var(iLookPARAM%zminLayer2)%dat(1),          & ! minimum layer depth for the 2nd layer (m)
 zminLayer3       => mpar_data%var(iLookPARAM%zminLayer3)%dat(1),          & ! minimum layer depth for the 3rd layer (m)
 zminLayer4       => mpar_data%var(iLookPARAM%zminLayer4)%dat(1),          & ! minimum layer depth for the 4th layer (m)
 zminLayer5       => mpar_data%var(iLookPARAM%zminLayer5)%dat(1),          & ! minimum layer depth for the 5th (bottom) layer (m)

 ! diagnostic scalar variables
 scalarSnowDepth  => prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),      & ! total snow depth (m)
 scalarSWE        => prog_data%var(iLookPROG%scalarSWE)%dat(1)             & ! SWE (kg m-2)

 ) ! end associate statement
 ! --------------------------------------------------------------------------------------------------------

 ! identify algorithmic control parameters to syb-divide and combine snow layers
 zminLayer = (/zminLayer1, zminLayer2, zminLayer3, zminLayer4, zminLayer5/)

 ! intialize the modified layers flag
 mergedLayers=.false.

 ! initialize the number of snow layers
 nSnow = indx_data%var(iLookINDEX%nSnow)%dat(1)
 nSoil   = indx_data%var(iLookINDEX%nSoil)%dat(1)
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)

 kSnow=0 ! initialize first layer to test (top layer)
 do ! attempt to remove multiple layers in a single time step (continuous do loop with exit clause)

  ! special case of >5 layers: add an offset to use maximum threshold from layer above
  if(ix_snowLayers == rulesDependLayerIndex .and. nSnow > 5)then
   nCheck=5
  else
   nCheck=nSnow
  end if

  ! loop through snow layers
  do iSnow=kSnow+1,nCheck

   ! associate local variables with the information in the data structures
   ! NOTE: do this here, since the layer variables are re-defined
   associate(&
   mLayerDepth      => prog_data%var(iLookPROG%mLayerDepth)%dat         , &    ! depth of each layer (m)
   mLayerVolFracIce => prog_data%var(iLookPROG%mLayerVolFracIce)%dat    , &    ! volumetric fraction of ice in each layer  (-)
   mLayerVolFracLiq => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat      &    ! volumetric fraction of liquid water in each layer (-)
   ) ! (associating local variables with the information in the data structures)

   ! check if the layer depth is less than the depth threshold
   select case(ix_snowLayers)
    case(sameRulesAllLayers);    removeLayer = (mLayerDepth(iSnow) < zmin)
    case(rulesDependLayerIndex); removeLayer = (mLayerDepth(iSnow) < zminLayer(iSnow))
    case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
   end select ! (option to combine/sub-divide snow layers)

   ! check if we have too much melt
   ! NOTE: assume that this is the top snow layer; need more trickery to relax this assumption
   if(tooMuchMelt .and. iSnow==1) removeLayer=.true.

   ! check if need to remove a layer
   if(removeLayer)then

    ! flag that we modified a layer
    mergedLayers=.true.

    ! ***** handle special case of a single layer
    if(nSnow==1)then
     ! set the variables defining "snow without a layer"
     ! NOTE: ignoring cold content!!! Need to fix later...
     scalarSnowDepth = mLayerDepth(1)
     scalarSWE       = (mLayerVolFracIce(1)*iden_ice + mLayerVolFracLiq(1)*iden_water)*mLayerDepth(1)
     ! remove the top layer from all model variable vectors
     ! NOTE: nSnow-1 = 0, so routine removes layer #1
     call rmLyAllVars(prog_data,prog_meta,nSnow-1,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     call rmLyAllVars(diag_data,diag_meta,nSnow-1,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     call rmLyAllVars(flux_data,flux_meta,nSnow-1,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     call rmLyAllVars(indx_data,indx_meta,nSnow-1,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if
     ! update the total number of layers
     nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
     nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)
     nLayers = nSnow + nSoil
     ! save the number of layers
     indx_data%var(iLookINDEX%nSnow)%dat(1)   = nSnow
     indx_data%var(iLookINDEX%nSoil)%dat(1)   = nSoil
     indx_data%var(iLookINDEX%nLayers)%dat(1) = nLayers
     ! update coordinate variables
     call calcHeight(&
                     ! input/output: data structures
                     indx_data,   & ! intent(in): layer type
                     prog_data,   & ! intent(inout): model variables for a local HRU
                     ! output: error control
                     err,cmessage)
     if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
     ! exit the do loop (no more snow layers to remove)
     return
    end if  ! (special case of 1 layer --> snow without a layer)

    ! ***** identify the layer to combine
    if(iSnow==1)then
     jSnow = iSnow+1  ! upper-most layer, combine with its lower neighbor
    elseif(iSnow==nSnow)then
     jSnow = nSnow-1  ! lower-most layer, combine with its upper neighbor
    else
     if(mLayerDepth(iSnow-1)<mLayerDepth(iSnow+1))then; jSnow = iSnow-1; else; jSnow = iSnow+1; end if
    end if

    ! ***** combine layers
    ! identify the layer closest to the surface
    kSnow=min(iSnow,jSnow)
    ! combine layer with identified neighbor
    call layer_combine(mpar_data,prog_data,diag_data,flux_data,indx_data,kSnow,err,cmessage)
    if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if

    ! exit the loop to try again
    exit

   end if  ! (if layer is below the mass threshold)

   kSnow=iSnow ! ksnow is used for completion test, so include here

   ! end association of local variables with the information in the data structures
   end associate

  end do ! (looping through snow layers)

  !print*, 'ksnow = ', ksnow

  ! exit if finished
  if(kSnow==nCheck)exit

 end do ! continuous do

 ! handle special case of > 5 layers in the CLM option
 if(nSnow > 5 .and. ix_snowLayers == rulesDependLayerIndex)then
  ! flag that layers were merged
  mergedLayers=.true.
  ! initial check to ensure everything is wonderful in the universe
  if(nSnow /= 6)then; err=5; message=trim(message)//'special case of > 5 layers: expect only six layers'; return; end if
  ! combine 5th layer with layer below
  call layer_combine(mpar_data,prog_data,diag_data,flux_data,indx_data,5,err,cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if
  ! another check
  if(nSnow /= 5)then; err=5; message=trim(message)//'special case of > 5 layers: expect to reduced layers to exactly 5'; return; end if
 end if

 ! check that there are no more than 5 layers in the CLM option
 if(ix_snowLayers == rulesDependLayerIndex)then
  if(nSnow > 5)then
   message=trim(message)//'expect no more than 5 layers when combination/sub-division rules depend on the layer index (CLM option)'
   err=20; return
  end if
 end if

 ! end association to variables in the data structure
 end associate

 end subroutine layerMerge


 ! ***********************************************************************************************************
 ! private subroutine layer_combine: combine snow layers and re-compute model state variables
 ! ***********************************************************************************************************
 ! combines layer iSnow with iSnow+1
 ! ***********************************************************************************************************
 subroutine layer_combine(mpar_data,prog_data,diag_data,flux_data,indx_data,iSnow,err,message)
 ! provide access to variables in the data structures
 USE var_lookup,only:iLookPARAM,iLookPROG,iLookINDEX           ! named variables for structure elements
 USE globalData,only:prog_meta,diag_meta,flux_meta,indx_meta   ! metadata
 USE data_types,only:var_ilength,var_dlength                   ! data vectors with variable length dimension
 USE data_types,only:var_d                                     ! data structures with fixed dimension
 ! provide access to external modules
 USE snow_utils_module,only:fracliquid                         ! compute fraction of liquid water
 USE convE2Temp_module,only:E2T_nosoil,temp2ethpy              ! convert temperature to enthalpy
 implicit none
 ! ------------------------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(var_dlength),intent(in)    :: mpar_data ! model parameters
 type(var_dlength),intent(inout) :: prog_data ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data ! model diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data ! model flux variables
 type(var_ilength),intent(inout) :: indx_data ! type of model layer
 ! input: snow layer indices
 integer(i4b),intent(in)         :: iSnow     ! index of top layer to combine
 ! output: error control
 integer(i4b),intent(out)        :: err       ! error code
 character(*),intent(out)        :: message   ! error message
 ! ------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=256)              :: cmessage                 ! error message for downwind routine
 real(dp)                        :: massIce(2)               ! mass of ice in the two layers identified for combination (kg m-2)
 real(dp)                        :: massLiq(2)               ! mass of liquid water in the two layers identified for combination (kg m-2)
 real(dp)                        :: bulkDenWat(2)            ! bulk density if total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(dp)                        :: cBulkDenWat              ! combined bulk density of total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(dp)                        :: cTemp                    ! combined layer temperature
 real(dp)                        :: cDepth                   ! combined layer depth
 real(dp)                        :: cVolFracIce              ! combined layer volumetric fraction of ice
 real(dp)                        :: cVolFracLiq              ! combined layer volumetric fraction of liquid water
 real(dp)                        :: l1Enthalpy,l2Enthalpy    ! enthalpy in the two layers identified for combination (J m-3)
 real(dp)                        :: cEnthalpy                ! combined layer enthalpy (J m-3)
 real(dp)                        :: fLiq                     ! fraction of liquid water at the combined temperature cTemp
 real(dp),parameter              :: eTol=1.e-4_dp            ! tolerance for the enthalpy-->temperature conversion (J m-3)

 ! initialize error control
 err=0; message="layer_combine/"

 ! associate local variables with information in the data structures
 associate(&
 ! model parameters
 snowfrz_scale    => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1), & ! scaling parameter for the freezing curve for snow (K-1)
 ! model state variables
 mLayerTemp       => prog_data%var(iLookPROG%mLayerTemp)%dat       , & ! temperature of each layer (K)
 mLayerDepth      => prog_data%var(iLookPROG%mLayerDepth)%dat      , & ! depth of each layer (m)
 mLayerVolFracIce => prog_data%var(iLookPROG%mLayerVolFracIce)%dat , & ! volumetric fraction of ice in each layer  (-)
 mLayerVolFracLiq => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat   & ! volumetric fraction of liquid water in each layer (-)
 ) ! (association of local variables with information in the data structures)

 ! compute combined depth
 cDepth       = mLayerDepth(isnow) + mLayerDepth(isnow+1)

 ! compute mass of each layer (kg m-2)
 massIce(1:2) = iden_ice*mLayerVolFracIce(iSnow:iSnow+1)*mLayerDepth(iSnow:iSnow+1)
 massLiq(1:2) = iden_water*mLayerVolFracLiq(iSnow:iSnow+1)*mLayerDepth(iSnow:iSnow+1)

 ! compute bulk density of water (kg m-3)
 bulkDenWat(1:2) = (massIce(1:2) + massLiq(1:2))/mLayerDepth(iSnow:iSnow+1)
 cBulkDenWat     = (mLayerDepth(isnow)*bulkDenWat(1) + mLayerDepth(isnow+1)*bulkDenWat(2))/cDepth

 ! compute enthalpy for each layer (J m-3)
 l1Enthalpy  = temp2ethpy(mLayerTemp(iSnow),  BulkDenWat(1),snowfrz_scale)
 l2Enthalpy  = temp2ethpy(mLayerTemp(iSnow+1),BulkDenWat(2),snowfrz_scale)

 ! compute combined enthalpy (J m-3)
 cEnthalpy   = (mLayerDepth(isnow)*l1Enthalpy + mLayerDepth(isnow+1)*l2Enthalpy)/cDepth

 ! convert enthalpy (J m-3) to temperature (K)
 call E2T_nosoil(cEnthalpy,cBulkDenWat,snowfrz_scale,cTemp,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if

 ! test enthalpy conversion
 if(abs(temp2ethpy(cTemp,cBulkDenWat,snowfrz_scale)/cBulkDenWat - cEnthalpy/cBulkDenWat) > eTol)then
  write(*,'(a,1x,f12.5,1x,2(e20.10,1x))') 'enthalpy test', cBulkDenWat, temp2ethpy(cTemp,cBulkDenWat,snowfrz_scale)/cBulkDenWat, cEnthalpy/cBulkDenWat
  message=trim(message)//'problem with enthalpy-->temperature conversion'
  err=20; return
 end if

 ! check temperature is within the two temperatures
 ! NOTE: use tolerance, for cases of merging a layer that has just been split
 if(cTemp > max(mLayerTemp(iSnow),mLayerTemp(iSnow+1))+eTol)then; err=20; message=trim(message)//'merged temperature > max(temp1,temp2)'; return; end if
 if(cTemp < min(mLayerTemp(iSnow),mLayerTemp(iSnow+1))-eTol)then; err=20; message=trim(message)//'merged temperature < min(temp1,temp2)'; return; end if

 ! compute volumetric fraction of liquid water
 fLiq = fracLiquid(cTemp,snowfrz_scale)

 ! compute volumetric fraction of ice and liquid water
 cVolFracLiq =          fLiq *cBulkDenWat/iden_water
 cVolFracIce = (1._dp - fLiq)*cBulkDenWat/iden_ice

 ! end association of local variables with information in the data structures
 end associate

 ! remove a model layer from all model variable vectors
 call rmLyAllVars(prog_data,prog_meta,iSnow,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 call rmLyAllVars(diag_data,diag_meta,iSnow,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 call rmLyAllVars(flux_data,flux_meta,iSnow,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 call rmLyAllVars(indx_data,indx_meta,iSnow,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! define the combined layer as snow
 indx_data%var(iLookINDEX%layerType)%dat(iSnow) = iname_snow

 ! update the total number of layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)
 nLayers = nSnow + nSoil

 ! save the number of layers in the data structures
 indx_data%var(iLookINDEX%nSnow)%dat(1)   = nSnow
 indx_data%var(iLookINDEX%nSoil)%dat(1)   = nSoil
 indx_data%var(iLookINDEX%nLayers)%dat(1) = nLayers

 ! ***** put state variables for the combined layer in the appropriate place
 prog_data%var(iLookPROG%mLayerTemp)%dat(iSnow)       = cTemp
 prog_data%var(iLookPROG%mLayerDepth)%dat(iSnow)      = cDepth
 prog_data%var(iLookPROG%mLayerVolFracIce)%dat(iSnow) = cVolFracIce
 prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(iSnow) = cVolFracLiq

 ! ***** adjust coordinate variables
 call calcHeight(&
                 ! input/output: data structures
                 indx_data,   & ! intent(in): layer type
                 prog_data,   & ! intent(inout): model variables for a local HRU
                 ! output: error control
                 err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 end subroutine layer_combine


 ! ***********************************************************************************************************
 ! private subroutine rmLyAllVars: reduce the length of the vectors in data structures
 ! ***********************************************************************************************************
 ! removes layer "iSnow+1" and sets layer "iSnow" to a missing value
 ! (layer "iSnow" will be filled with a combined layer later)
 ! ***********************************************************************************************************
 subroutine rmLyAllVars(dataStruct,metaStruct,iSnow,err,message)
 USE var_lookup,only:iLookVarType                 ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName       ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc            ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength      ! data vectors with variable length dimension
 USE data_types,only:var_info                     ! metadata structure
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 class(*),intent(inout)          :: dataStruct     ! data structure
 type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in)         :: iSnow          ! new layer
 ! output: error control
 integer(i4b),intent(out)        :: err            ! error code
 character(*),intent(out)        :: message        ! error message
 ! locals
 integer(i4b)                    :: ivar           ! variable index
 integer(i4b)                    :: ix_lower       ! lower bound of the vector
 integer(i4b)                    :: ix_upper       ! upper bound of the vector
 real(dp),allocatable            :: tempVec_dp(:)  ! temporary vector (double precision)
 integer(i4b),allocatable        :: tempVec_i4b(:) ! temporary vector (integer)
 character(LEN=256)              :: cmessage       ! error message of downwind routine
 ! initialize error control
 err=0; message="rmLyAllVars/"

 ! check dimensions
 select type(dataStruct)
  type is (var_dlength); if(size(dataStruct%var) /= size(metaStruct)) err=20
  type is (var_ilength); if(size(dataStruct%var) /= size(metaStruct)) err=20
  class default; err=20; message=trim(message)//'unable to identify the data type'; return
 end select
 if(err/=0)then; message=trim(message)//'dimensions of data structure and metadata structures do not match'; return; end if

 ! ***** loop through model variables and remove one layer
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=nSnow
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=nSnow
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nLayers
   case default; cycle  ! no need to remove soil layers or scalar variables
  end select

  ! remove layers
  select type(dataStruct)

   ! ** double precision
   type is (var_dlength)
    ! check allocated
    if(.not.allocated(dataStruct%var(ivar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! allocate the temporary vector
    allocate(tempVec_dp(ix_lower:ix_upper-1), stat=err)
    if(err/=0)then; err=20; message=trim(message)//'unable to allocate temporary vector'; return; end if
    ! copy elements across to the temporary vector
    if(iSnow>=ix_lower)  tempVec_dp(iSnow)              = missingDouble ! set merged layer to missing (fill in later)
    if(iSnow>ix_lower)   tempVec_dp(ix_lower:iSnow-1)   = dataStruct%var(ivar)%dat(ix_lower:iSnow-1)
    if(iSnow+1<ix_upper) tempVec_dp(iSnow+1:ix_upper-1) = dataStruct%var(ivar)%dat(iSnow+2:ix_upper)  ! skip iSnow+1
    ! deallocate the data vector: strictly not necessary, but include to be safe 
    deallocate(dataStruct%var(ivar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem deallocating data vector'; return; end if
    ! create the new data structure using the temporary vector as the source
    call cloneStruc(dataStruct%var(ivar)%dat, ix_lower, source=tempVec_dp, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! deallocate the temporary data vector: strictly not necessary, but include to be safe
    deallocate(tempVec_dp,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! ** integer
   type is (var_ilength)
    ! check allocated
    if(.not.allocated(dataStruct%var(ivar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! allocate the temporary vector
    allocate(tempVec_i4b(ix_lower:ix_upper-1), stat=err)
    if(err/=0)then; err=20; message=trim(message)//'unable to allocate temporary vector'; return; end if
    ! copy elements across to the temporary vector
    if(iSnow>=ix_lower)  tempVec_i4b(iSnow)              = missingInteger ! set merged layer to missing (fill in later)
    if(iSnow>ix_lower)   tempVec_i4b(ix_lower:iSnow-1)   = dataStruct%var(ivar)%dat(ix_lower:iSnow-1)
    if(iSnow+1<ix_upper) tempVec_i4b(iSnow+1:ix_upper-1) = dataStruct%var(ivar)%dat(iSnow+2:ix_upper)  ! skip iSnow+1
    ! deallocate the data vector: strictly not necessary, but include to be safe
    deallocate(dataStruct%var(ivar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem deallocating data vector'; return; end if
    ! create the new data structure using the temporary vector as the source
    call cloneStruc(dataStruct%var(ivar)%dat, ix_lower, source=tempVec_i4b, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! deallocate the temporary data vector: strictly not necessary, but include to be safe
    deallocate(tempVec_i4b,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! check that we found the data type
   class default; err=20; message=trim(message)//'unable to identify the data type'; return

  end select ! dependence on data types

 end do  ! looping through variables

 end subroutine rmLyAllVars

end module layerMerge_module
