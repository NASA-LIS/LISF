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

module ssdNrgFlux_module
! numerical recipes data types
USE nrtype
! physical constants
USE multiconst,only:&
                    sb,          & ! Stefan Boltzman constant      (W m-2 K-4)
                    Em_Sno,      & ! emissivity of snow            (-)
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_water,    & ! specifric heat of water       (J kg-1 K-1)
                    LH_fus,      & ! latent heat of fusion         (J kg-1)
                    LH_vap,      & ! latent heat of vaporization   (J kg-1)
                    LH_sub,      & ! latent heat of sublimation    (J kg-1)
                    gravity,     & ! gravitational acceleteration  (m s-2)
                    Tfreeze,     & ! freezing point of pure water  (K)
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water     ! intrinsic density of water    (kg m-3)
! named variables for snow and soil
USE globalData,only:iname_snow     ! named variables for snow
USE globalData,only:iname_soil     ! named variables for soil
! provide access to look-up values for model decisions
USE mDecisions_module,only:      &
 ! look-up values for the numerical method
 iterative,                      & ! iterative
 nonIterative,                   & ! non-iterative
 iterSurfEnergyBal,              & ! iterate only on the surface energy balance
 ! look-up values for method used to compute derivative
 numerical,                      & ! numerical solution
 analytical,                     & ! analytical solution
 ! look-up values for choice of boundary conditions for thermodynamics
 prescribedTemp,                 & ! prescribed temperature
 energyFlux,                     & ! energy flux
 zeroFlux,                       & ! zero flux
 ! look-up values for choice of boundary conditions for soil hydrology
 prescribedHead                    ! prescribed head
! -------------------------------------------------------------------------------------------------
implicit none
private
public::ssdNrgFlux
! global parameters
real(dp),parameter            :: dx=1.e-10_dp             ! finite difference increment (K)
real(dp),parameter            :: valueMissing=-9999._dp   ! missing value parameter
contains


 ! ************************************************************************************************
 ! public subroutine ssdNrgFlux: compute energy fluxes and derivatives at layer interfaces
 ! ************************************************************************************************
 subroutine ssdNrgFlux(&
                       ! input: fluxes and derivatives at the upper boundary
                       groundNetFlux,                      & ! intent(in):    total flux at the ground surface (W m-2)
                       dGroundNetFlux_dGroundTemp,         & ! intent(in):    derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                       ! input: liquid water fluxes
                       iLayerLiqFluxSnow,                  & ! intent(in):    liquid flux at the interface of each snow layer (m s-1)
                       iLayerLiqFluxSoil,                  & ! intent(in):    liquid flux at the interface of each soil layer (m s-1)
                       ! input: trial value of model state variabes
                       mLayerTempTrial,                    & ! intent(in):    trial temperature at the current iteration (K)
                       ! input-output: data structures
                       mpar_data,                          & ! intent(in):    model parameters
                       indx_data,                          & ! intent(in):    model indices
                       prog_data,                          & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,                          & ! intent(in):    model diagnostic variables for a local HRU
                       flux_data,                          & ! intent(inout): model fluxes for a local HRU
                       ! output: fluxes and derivatives at all layer interfaces
                       iLayerNrgFlux,                      & ! intent(out):   energy flux at the layer interfaces (W m-2)
                       dFlux_dTempAbove,                   & ! intent(out):   derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                       dFlux_dTempBelow,                   & ! intent(out):   derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                       ! output: error control
                       err,message)                          ! intent(out): error control
 ! model decisions
 USE globalData,only:model_decisions                         ! model decision structure
 USE var_lookup,only:iLookDECISIONS                          ! named variables for elements of the decision structure
 ! named variables
 USE var_lookup,only:iLookPROG       ! named variables for structure elements
 USE var_lookup,only:iLookDIAG       ! named variables for structure elements
 USE var_lookup,only:iLookFLUX       ! named variables for structure elements
 USE var_lookup,only:iLookPARAM      ! named variables for structure elements
 USE var_lookup,only:iLookINDEX      ! named variables for structure elements
 ! data types
 USE data_types,only:var_d           ! x%var(:)       (dp)
 USE data_types,only:var_ilength     ! x%var(:)%dat   (i4b)
 USE data_types,only:var_dlength     ! x%var(:)%dat   (dp)
 implicit none
 ! input: fluxes and derivatives at the upper boundary
 real(dp),intent(in)             :: groundNetFlux              ! net energy flux for the ground surface (W m-2)
 real(dp),intent(in)             :: dGroundNetFlux_dGroundTemp ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
 ! input: liquid water fluxes
 real(dp),intent(in)             :: iLayerLiqFluxSnow(0:)      ! intent(in): liquid flux at the interface of each snow layer (m s-1)
 real(dp),intent(in)             :: iLayerLiqFluxSoil(0:)      ! intent(in): liquid flux at the interface of each soil layer (m s-1)
 ! input: trial value of model state variables
 real(dp),intent(in)             :: mLayerTempTrial(:)         ! trial temperature of each snow/soil layer at the current iteration (K)
 ! input-output: data structures
 type(var_dlength),intent(in)    :: mpar_data                  ! model parameters
 type(var_ilength),intent(in)    :: indx_data                  ! state vector geometry
 type(var_dlength),intent(in)    :: prog_data                  ! prognostic variables for a local HRU
 type(var_dlength),intent(in)    :: diag_data                  ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                  ! model fluxes for a local HRU
 ! output: fluxes and derivatives at all layer interfaces
 real(dp),intent(out)            :: iLayerNrgFlux(0:)          ! energy flux at the layer interfaces (W m-2)
 real(dp),intent(out)            :: dFlux_dTempAbove(0:)       ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),intent(out)            :: dFlux_dTempBelow(0:)       ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 ! output: error control
 integer(i4b),intent(out)        :: err                        ! error code
 character(*),intent(out)        :: message                    ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: iLayer                     ! index of model layers
 real(dp)                        :: qFlux                      ! liquid flux at layer interfaces (m s-1)
 real(dp)                        :: dz                         ! height difference (m)
 real(dp)                        :: flux0,flux1,flux2          ! fluxes used to calculate derivatives (W m-2)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! make association of local variables with information in the data structures
 associate(&
  ix_fDerivMeth        => model_decisions(iLookDECISIONS%fDerivMeth)%iDecision, & ! intent(in): method used to calculate flux derivatives
  ix_bcLowrTdyn        => model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision, & ! intent(in): method used to calculate the lower boundary condition for thermodynamics
  ! input: model coordinates and thermal properties
  nSnow                => indx_data%var(iLookINDEX%nSnow)%dat(1),               & ! intent(in): number of snow layers 
  nLayers              => indx_data%var(iLookINDEX%nLayers)%dat(1),             & ! intent(in): total number of layers 
  layerType            => indx_data%var(iLookINDEX%layerType)%dat,              & ! intent(in): layer type (iname_soil or iname_snow)
  mLayerDepth          => prog_data%var(iLookPROG%mLayerDepth)%dat,             & ! intent(in): depth of each layer (m)
  mLayerHeight         => prog_data%var(iLookPROG%mLayerHeight)%dat,            & ! intent(in): height at the mid-point of each layer (m)
  iLayerThermalC       => diag_data%var(iLookDIAG%iLayerThermalC)%dat,          & ! intent(in): thermal conductivity at the interface of each layer (W m-1 K-1)
  lowerBoundTemp       => mpar_data%var(iLookPARAM%lowerBoundTemp)%dat(1),      & ! intent(in): temperature of the lower boundary (K)
  ! output: diagnostic fluxes
  iLayerConductiveFlux => flux_data%var(iLookFLUX%iLayerConductiveFlux)%dat,    & ! intent(out): conductive energy flux at layer interfaces at end of time step (W m-2)
  iLayerAdvectiveFlux  => flux_data%var(iLookFLUX%iLayerAdvectiveFlux)%dat      & ! intent(out): advective energy flux at layer interfaces at end of time step (W m-2)
 )  ! association of local variables with information in the data structures 
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='ssdNrgFlux/'

 ! set conductive and advective fluxes to missing in the upper boundary
 ! NOTE: advective flux at the upper boundary is included in the ground heat flux
 iLayerConductiveFlux(0) = valueMissing
 iLayerAdvectiveFlux(0)  = valueMissing

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the conductive fluxes at layer interfaces *****
 ! -------------------------------------------------------------------------------------------------------------------------
 do iLayer=1,nLayers

  ! compute fluxes at the lower boundary -- positive downwards
  if(iLayer==nLayers)then
   ! flux depends on the type of lower boundary condition
   select case(ix_bcLowrTdyn) ! (identify the lower boundary condition for thermodynamics
    case(prescribedTemp); iLayerConductiveFlux(nLayers) = -iLayerThermalC(iLayer)*(lowerBoundTemp - mLayerTempTrial(iLayer))/(mLayerDepth(iLayer)*0.5_dp)
    case(zeroFlux);       iLayerConductiveFlux(nLayers) = 0._dp
    case default;         err=20; message=trim(message)//'unable to identify lower boundary condition for thermodynamics'; return
   end select  ! (identifying the lower boundary condition for thermodynamics)

  ! compute fluxes within the domain -- positive downwards
  else
    iLayerConductiveFlux(iLayer)  = -iLayerThermalC(iLayer)*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer)) / &
                                    (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))

    !write(*,'(a,i4,1x,2(f9.3,1x))') 'iLayer, iLayerConductiveFlux(iLayer), iLayerThermalC(iLayer) = ', iLayer, iLayerConductiveFlux(iLayer), iLayerThermalC(iLayer)
  end if ! (the type of layer)
 end do

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the advective fluxes at layer interfaces *****
 ! -------------------------------------------------------------------------------------------------------------------------
 do iLayer=1,nLayers
  ! get the liquid flux at layer interfaces
  select case(layerType(iLayer))
   case(iname_snow); qFlux = iLayerLiqFluxSnow(iLayer)
   case(iname_soil); qFlux = iLayerLiqFluxSoil(iLayer-nSnow)
   case default; err=20; message=trim(message)//'unable to identify layer type'; return
  end select
  ! compute fluxes at the lower boundary -- positive downwards
  if(iLayer==nLayers)then
   iLayerAdvectiveFlux(iLayer) = -Cp_water*iden_water*qFlux*(lowerBoundTemp - mLayerTempTrial(iLayer))
  ! compute fluxes within the domain -- positive downwards
  else
   iLayerAdvectiveFlux(iLayer) = -Cp_water*iden_water*qFlux*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer))
  end if
 end do  ! looping through layers

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the total fluxes at layer interfaces *****
 ! -------------------------------------------------------------------------------------------------------------------------
 ! NOTE: ignore advective fluxes for now
 iLayerNrgFlux(0)         = groundNetFlux
 iLayerNrgFlux(1:nLayers) = iLayerConductiveFlux(1:nLayers)
 !print*, 'iLayerNrgFlux(0:4) = ', iLayerNrgFlux(0:4)

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the derivative in fluxes at layer interfaces w.r.t temperature in the layer above and the layer below *****
 ! -------------------------------------------------------------------------------------------------------------------------

 ! initialize un-used elements
 dFlux_dTempBelow(nLayers) = -huge(lowerBoundTemp)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems

 ! loop through INTERFACES...
 do iLayer=0,nLayers

  ! ***** the upper boundary -- ** NOTE: dTotalSurfaceFlux_dTemp was computed previously using ix_fDerivMeth
  if(iLayer==0)then  ! (upper boundary)
   dFlux_dTempBelow(iLayer) = dGroundNetFlux_dGroundTemp

  ! ***** the lower boundary
  elseif(iLayer==nLayers)then  ! (lower boundary)

   ! identify the lower boundary condition
   select case(ix_bcLowrTdyn)

    ! * prescribed temperature at the lower boundary
    case(prescribedTemp)

     dz = mLayerDepth(iLayer)*0.5_dp
     if(ix_fDerivMeth==analytical)then    ! ** analytical derivatives
      dFlux_dTempAbove(iLayer) = iLayerThermalC(iLayer)/dz
     else                              ! ** numerical derivatives
      flux0 = -iLayerThermalC(iLayer)*(lowerBoundTemp - (mLayerTempTrial(iLayer)   ))/dz
      flux1 = -iLayerThermalC(iLayer)*(lowerBoundTemp - (mLayerTempTrial(iLayer)+dx))/dz
      dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
     end if

     ! * zero flux at the lower boundary
     case(zeroFlux)
      dFlux_dTempAbove(iLayer) = 0._dp

     case default; err=20; message=trim(message)//'unable to identify lower boundary condition for thermodynamics'; return

   end select  ! (identifying the lower boundary condition for thermodynamics)

  ! ***** internal layers
  else
   dz = (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
   if(ix_fDerivMeth==analytical)then    ! ** analytical derivatives
    dFlux_dTempAbove(iLayer) =  iLayerThermalC(iLayer)/dz
    dFlux_dTempBelow(iLayer) = -iLayerThermalC(iLayer)/dz
   else                              ! ** numerical derivatives
    flux0 = -iLayerThermalC(iLayer)*( mLayerTempTrial(iLayer+1)     -  mLayerTempTrial(iLayer)    ) / dz
    flux1 = -iLayerThermalC(iLayer)*( mLayerTempTrial(iLayer+1)     - (mLayerTempTrial(iLayer)+dx)) / dz
    flux2 = -iLayerThermalC(iLayer)*((mLayerTempTrial(iLayer+1)+dx) -  mLayerTempTrial(iLayer)    ) / dz
    dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
    dFlux_dTempBelow(iLayer) = (flux2 - flux0)/dx
   end if

  end if  ! type of layer (upper, internal, or lower)

 end do  ! (looping through layers)

 ! end association of local variables with information in the data structures
 end associate

 end subroutine ssdNrgFlux

end module ssdNrgFlux_module

