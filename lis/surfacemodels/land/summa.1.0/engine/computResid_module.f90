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

module computResid_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! define access to state variables to print
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (dp)

implicit none
private
public::computResid
contains

 ! **********************************************************************************************************
 ! public subroutine computResid: compute the residual vector
 ! **********************************************************************************************************
 subroutine computResid(&
                        ! input: model control
                        dt,                      & ! intent(in):    length of the time step (seconds)
                        nSnow,                   & ! intent(in):    number of snow layers
                        nSoil,                   & ! intent(in):    number of soil layers
                        nLayers,                 & ! intent(in):    total number of layers
                        ! input: flux vectors
                        sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                        fVec,                    & ! intent(in):    flux vector
                        ! input: state variables (already disaggregated into scalars and vectors)
                        scalarCanairTempTrial,   & ! intent(in):    trial value for the temperature of the canopy air space (K)
                        scalarCanopyTempTrial,   & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                        scalarCanopyHydTrial,    & ! intent(in):    trial value of canopy water (kg m-2), either liquid water content or total water content
                        mLayerTempTrial,         & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                        mLayerVolFracHydTrial,   & ! intent(in):    trial vector of volumetric water content (-), either liquid water content or total water content
                        ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                        scalarCanopyIceTrial,    & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                        mLayerVolFracIceTrial,   & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                        ! input: data structures
                        prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                        flux_data,               & ! intent(in):    model fluxes for a local HRU
                        indx_data,               & ! intent(in):    index data
                        ! output
                        rAdd,                    & ! intent(out):   additional (sink) terms on the RHS of the state equation
                        rVec,                    & ! intent(out):   residual vector
                        err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 USE var_lookup,only:iLookPROG                     ! named variables for structure elements
 USE var_lookup,only:iLookDIAG                     ! named variables for structure elements
 USE var_lookup,only:iLookFLUX                     ! named variables for structure elements
 USE var_lookup,only:iLookINDEX                    ! named variables for structure elements
 implicit none
 ! input: model control
 real(dp),intent(in)             :: dt                        ! length of the time step (seconds)
 integer(i4b),intent(in)         :: nSnow                     ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                     ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                   ! total number of layers in the snow+soil domain
 ! input: flux vectors
 real(qp),intent(in)             :: sMul(:)   ! NOTE: qp      ! state vector multiplier (used in the residual calculations)
 real(dp),intent(in)             :: fVec(:)                   ! flux vector
 ! input: state variables (already disaggregated into scalars and vectors)
 real(dp),intent(in)             :: scalarCanairTempTrial     ! trial value for temperature of the canopy air space (K)
 real(dp),intent(in)             :: scalarCanopyTempTrial     ! trial value for temperature of the vegetation canopy (K)
 real(dp),intent(in)             :: scalarCanopyHydTrial      ! trial value for canopy water (kg m-2), either liquid water content or total water content
 real(dp),intent(in)             :: mLayerTempTrial(:)        ! trial value for temperature of each snow/soil layer (K)
 real(dp),intent(in)             :: mLayerVolFracHydTrial(:)  ! trial vector of volumetric water content (-), either liquid water content or total water content 
 ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
 real(dp),intent(in)             :: scalarCanopyIceTrial      ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)             :: mLayerVolFracIceTrial(:)  ! trial value for volumetric fraction of ice (-)
 ! input: data structures
 type(var_dlength),intent(in)    :: prog_data                 ! prognostic variables for a local HRU
 type(var_dlength),intent(in)    :: diag_data                 ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)    :: flux_data                 ! model fluxes for a local HRU
 type(var_ilength),intent(in)    :: indx_data                 ! indices defining model states and layers
 ! output
 real(dp),intent(out)            :: rAdd(:)                   ! additional (sink) terms on the RHS of the state equation
 real(qp),intent(out)            :: rVec(:)   ! NOTE: qp      ! residual vector
 integer(i4b),intent(out)        :: err                       ! error code
 character(*),intent(out)        :: message                   ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 integer(i4b)                    :: iLayer                    ! index of layer within the snow+soil domain
 integer(i4b),parameter          :: ixVegVolume=1             ! index of the desired vegetation control volumne (currently only one veg layer)
 real(dp)                        :: scalarCanopyHyd           ! canopy water content (kg m-2), either liquid water content or total water content
 real(dp),dimension(nLayers)     :: mLayerVolFracHyd          ! vector of volumetric water content (-), either liquid water content or total water content
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! link to the necessary variables for the residual computations
 associate(&
  ! model state variables (vegetation canopy)
  scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in): [dp]     temperature of the canopy air space (K)
  scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in): [dp]     temperature of the vegetation canopy (K)
  scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(in): [dp]     mass of ice on the vegetation canopy (kg m-2)
  scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(in): [dp]     mass of liquid water on the vegetation canopy (kg m-2)
  scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in): [dp]     mass of total water on the vegetation canopy (kg m-2)
  ! model state variables (snow and soil domains)
  mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in): [dp(:)]  temperature of each snow/soil layer (K)
  mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(in): [dp(:)]  volumetric fraction of ice (-)
  mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in): [dp(:)]  volumetric fraction of liquid water (-)
  mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in): [dp(:)]  volumetric fraction of total water (-)
  ! canopy and layer depth
  canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in): [dp]      canopy depth (m)
  mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
  ! model fluxes (sink terms in the soil domain)
  mLayerTranspire         => flux_data%var(iLookFLUX%mLayerTranspire)%dat           ,& ! intent(in): [dp]     transpiration loss from each soil layer (m s-1)
  mLayerBaseflow          => flux_data%var(iLookFLUX%mLayerBaseflow)%dat            ,& ! intent(in): [dp(:)]  baseflow from each soil layer (m s-1)
  mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,& ! intent(in): [dp(:)]  change in storage associated with compression of the soil matrix (-)
  ! number of state variables of a specific type
  nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
  nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
  nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
  ! model indices
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy energy state variable
  ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
  ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
  ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
  ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
  ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in): [i4b(:)] named variables defining the type of hydrology states in snow+soil domain
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                 & ! intent(in): [i4b(:)] named variables defining the type of layer in snow+soil domain
 ) ! association to necessary variables for the residual computations
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="computResid/"

 ! ---
 ! * compute sink terms...
 ! -----------------------

 ! intialize additional terms on the RHS as zero
 rAdd(:) = 0._dp

 ! compute energy associated with melt freeze for the vegetation canopy
 if(ixVegNrg/=integerMissing) rAdd(ixVegNrg) = rAdd(ixVegNrg) + LH_fus*(scalarCanopyIceTrial - scalarCanopyIce)/canopyDepth   ! energy associated with melt/freeze (J m-3)

 ! compute energy associated with melt/freeze for snow
 ! NOTE: allow expansion of ice during melt-freeze for snow; deny expansion of ice during melt-freeze for soil
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   select case( layerType(iLayer) )
    case(iname_snow); rAdd( ixSnowSoilNrg(iLayer) ) = rAdd( ixSnowSoilNrg(iLayer) ) + LH_fus*iden_ice  *(mLayerVolFracIceTrial(iLayer) - mLayerVolFracIce(iLayer))
    case(iname_soil); rAdd( ixSnowSoilNrg(iLayer) ) = rAdd( ixSnowSoilNrg(iLayer) ) + LH_fus*iden_water*(mLayerVolFracIceTrial(iLayer) - mLayerVolFracIce(iLayer))
   end select
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! sink terms soil hydrology (-)
 ! NOTE 1: state variable is volumetric water content, so melt-freeze is not included
 ! NOTE 2: ground evaporation was already included in the flux at the upper boundary
 ! NOTE 3: rAdd(ixSnowOnlyWat)=0, and is defined in the initialization above
 ! NOTE 4: same sink terms for matric head and liquid matric potential
 if(nSoilOnlyHyd>0)then
  do concurrent (iLayer=1:nSoil,ixSoilOnlyHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
   rAdd( ixSoilOnlyHyd(iLayer) ) = rAdd( ixSoilOnlyHyd(iLayer) ) + dt*(mLayerTranspire(iLayer) - mLayerBaseflow(iLayer) )/mLayerDepth(iLayer+nSnow) - mLayerCompress(iLayer)
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! ---
 ! * compute the residual vector...
 ! --------------------------------

 ! compute the residual vector for the vegetation canopy
 ! NOTE: sMul(ixVegHyd) = 1, but include as it converts all variables to quadruple precision
 ! --> energy balance
 if(ixCasNrg/=integerMissing) rVec(ixCasNrg) = sMul(ixCasNrg)*scalarCanairTempTrial - ( (sMul(ixCasNrg)*scalarCanairTemp + fVec(ixCasNrg)*dt) + rAdd(ixCasNrg) )
 if(ixVegNrg/=integerMissing) rVec(ixVegNrg) = sMul(ixVegNrg)*scalarCanopyTempTrial - ( (sMul(ixVegNrg)*scalarCanopyTemp + fVec(ixVegNrg)*dt) + rAdd(ixVegNrg) )
 ! --> mass balance
 if(ixVegHyd/=integerMissing)then
  scalarCanopyHyd = merge(scalarCanopyWat, scalarCanopyLiq, (ixStateType( ixHydCanopy(ixVegVolume) )==iname_watCanopy) )
  rVec(ixVegHyd)  = sMul(ixVegHyd)*scalarCanopyHydTrial  - ( (sMul(ixVegHyd)*scalarCanopyHyd  + fVec(ixVegHyd)*dt) + rAdd(ixVegHyd) )
 endif

 ! compute the residual vector for the snow and soil sub-domains for energy
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   rVec( ixSnowSoilNrg(iLayer) ) = sMul( ixSnowSoilNrg(iLayer) )*mLayerTempTrial(iLayer) - ( (sMul( ixSnowSoilNrg(iLayer) )*mLayerTemp(iLayer)  + fVec( ixSnowSoilNrg(iLayer) )*dt) + rAdd( ixSnowSoilNrg(iLayer) ) )
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! compute the residual vector for the snow and soil sub-domains for hydrology
 ! NOTE: residual depends on choice of state variable
 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
   ! (get the correct state variable)
   mLayerVolFracHyd(iLayer)      = merge(mLayerVolFracWat(iLayer), mLayerVolFracLiq(iLayer), (ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_matLayer) )
   ! (compute the residual)
   rVec( ixSnowSoilHyd(iLayer) ) = mLayerVolFracHydTrial(iLayer) - ( (mLayerVolFracHyd(iLayer) + fVec( ixSnowSoilHyd(iLayer) )*dt) + rAdd( ixSnowSoilHyd(iLayer) ) )
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! print result
 if(globalPrintFlag)then
  write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
  write(*,'(a,1x,100(e12.5,1x))') 'fVec = ', fVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
  !print*, 'PAUSE:'; read(*,*)
 endif

 ! check
 if(any(isNan(rVec)))then
  message=trim(message)//'we found some Indian bread (NaN)'
  err=20; return
 endif

 ! end association with the necessary variabiles for the residual calculations
 end associate

 end subroutine computResid

end module computResid_module
