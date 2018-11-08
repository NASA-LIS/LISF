module flxMapping_module
implicit none
private
public::flxMapping
contains

 subroutine flxMapping(err,message)
 USE nrtype
 ! data types
 USE data_types, only: var_info        ! data type for metadata structure
 USE data_types, only: flux2state      ! data type for extended metadata structure, for flux-to-state mapping
 ! structures of named variables
 USE var_lookup, only: iLookFLUX       ! named variables for local flux variables
 ! metadata structures
 USE globalData, only: flux_meta       ! data structure for model fluxes
 USE globalData, only: flux2state_orig ! data structure for flux-to-state mapping (original state variables)
 USE globalData, only: flux2state_liq  ! data structure for flux-to-state mapping (liquid water state variables)
 ! named variables to describe the state variable type
 USE globalData, only: iname_nrgCanair ! named variable defining the energy of the canopy air space
 USE globalData, only: iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
 USE globalData, only: iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
 USE globalData, only: iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
 USE globalData, only: iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
 USE globalData, only: iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
 USE globalData, only: iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
 USE globalData, only: iname_matLayer  ! named variable defining the matric head state variable for soil layers
 USE globalData, only: iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers
 ! access missing values
 USE globalData,only:integerMissing    ! missing integer
 implicit none
 ! dummy variables
 integer(i4b),intent(out)       :: err                 ! error code
 character(*),intent(out)       :: message             ! error message
 ! local variables
 integer(i4b)                   :: iVar                ! variable index
 integer(i4b)                   :: nFlux               ! number of fluxes
 integer(i4b),parameter         :: integerUndefined=0  ! named variable to denote that the flux is undefined
 ! initialize error control
 err=0; message='flxMapping/'

 ! get the number of fluxes
 nFlux = size(flux_meta)

 ! -----
 ! - original state variables...
 ! -----------------------------

 ! ** initialize flux-to-state mapping
 do iVar=1,nFlux
  flux2state_orig(iVar)%state1 = integerUndefined
  flux2state_orig(iVar)%state2 = integerUndefined
 end do

 ! ** define mapping between fluxes and states

 ! net energy and mass fluxes for the vegetation domain
 flux2state_orig(iLookFLUX%scalarCanopyNetLiqFlux)          = flux2state(state1=iname_watCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarCanairNetNrgFlux)          = flux2state(state1=iname_nrgCanair, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarCanopyNetNrgFlux)          = flux2state(state1=iname_nrgCanopy, state2=integerMissing) 
 flux2state_orig(iLookFLUX%scalarGroundNetNrgFlux)          = flux2state(state1=iname_nrgLayer,  state2=integerMissing)

 ! precipitation -- does not depend on state variables
 flux2state_orig(iLookFLUX%scalarRainfall)                  = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarSnowfall)                  = flux2state(state1=integerMissing, state2=integerMissing)
 
 ! shortwave radiation -- does not depend on state variables
 flux2state_orig(iLookFLUX%spectralIncomingDirect)          = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_orig(iLookFLUX%spectralIncomingDiffuse)         = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarCanopySunlitPAR)           = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarCanopyShadedPAR)           = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_orig(iLookFLUX%spectralBelowCanopyDirect)       = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_orig(iLookFLUX%spectralBelowCanopyDiffuse)      = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarBelowCanopySolar)          = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarCanopyAbsorbedSolar)       = flux2state(state1=integerMissing, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarGroundAbsorbedSolar)       = flux2state(state1=integerMissing, state2=integerMissing)
 
 ! longwave radiation -- assume calculated when the canopy energy state variable is active OR when the ground energy state variable is active
 flux2state_orig(iLookFLUX%scalarLWRadCanopy)               = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarLWRadGround)               = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarLWRadUbound2Canopy)        = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarLWRadUbound2Ground)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarLWRadUbound2Ubound)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarLWRadCanopy2Ubound)        = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarLWRadCanopy2Ground)        = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarLWRadCanopy2Canopy)        = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarLWRadGround2Ubound)        = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarLWRadGround2Canopy)        = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarLWNetCanopy)               = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarLWNetGround)               = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarLWNetUbound)               = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 
 ! turbulent heat transfer -- assume calculated when the canopy energy state variable is active OR when the ground energy state variable is active 
 flux2state_orig(iLookFLUX%scalarEddyDiffusCanopyTop)       = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarFrictionVelocity)          = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarWindspdCanopyTop)          = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarWindspdCanopyBottom)       = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarGroundResistance)          = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarCanopyResistance)          = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarLeafResistance)            = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarSoilResistance)            = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarSenHeatTotal)              = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarSenHeatCanopy)             = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarSenHeatGround)             = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarLatHeatTotal)              = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarLatHeatCanopyEvap)         = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarLatHeatCanopyTrans)        = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarLatHeatGround)             = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarCanopyAdvectiveHeatFlux)   = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarGroundAdvectiveHeatFlux)   = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarCanopySublimation)         = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarSnowSublimation)           = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 
 ! stomatal resistance and photosynthesis -- calculated when the canopy energy state variable is active
 flux2state_orig(iLookFLUX%scalarStomResistSunlit)          = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarStomResistShaded)          = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarPhotosynthesisSunlit)      = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarPhotosynthesisShaded)      = flux2state(state1=iname_nrgCanopy, state2=integerMissing)

 ! liquid water fluxes associated with evapotranspiration
 ! NOTE 1: calculated in the energy balance routines: energy balance must be calculated first in order for water to balance
 ! NOTE 2: if implement strang splitting, need to average fluxes from the start and end of the time step
 flux2state_orig(iLookFLUX%scalarCanopyTranspiration)       = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%scalarCanopyEvaporation)         = flux2state(state1=iname_nrgCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarGroundEvaporation)         = flux2state(state1=iname_nrgCanopy, state2=iname_nrgLayer)
 flux2state_orig(iLookFLUX%mLayerTranspire)                 = flux2state(state1=iname_matLayer,  state2=integerMissing)

 ! liquid and solid water fluxes through the canopy
 flux2state_orig(iLookFLUX%scalarThroughfallSnow)           = flux2state(state1=integerMissing,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarCanopySnowUnloading)       = flux2state(state1=integerMissing,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarThroughfallRain)           = flux2state(state1=iname_watCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarCanopyLiqDrainage)         = flux2state(state1=iname_watCanopy, state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarCanopyMeltFreeze)          = flux2state(state1=integerMissing,  state2=integerMissing)
 
 ! energy fluxes and for the snow and soil domains
 flux2state_orig(iLookFLUX%iLayerConductiveFlux)            = flux2state(state1=iname_nrgLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%iLayerAdvectiveFlux)             = flux2state(state1=iname_nrgLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%iLayerNrgFlux)                   = flux2state(state1=iname_nrgLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%mLayerNrgFlux)                   = flux2state(state1=iname_nrgLayer,  state2=integerMissing)
 
 ! liquid water fluxes for the snow domain
 flux2state_orig(iLookFLUX%scalarSnowDrainage)              = flux2state(state1=iname_watLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%iLayerLiqFluxSnow)               = flux2state(state1=iname_watLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%mLayerLiqFluxSnow)               = flux2state(state1=iname_watLayer,  state2=integerMissing)
 
 ! liquid water fluxes for the soil domain
 flux2state_orig(iLookFLUX%scalarRainPlusMelt)              = flux2state(state1=iname_watLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarMaxInfilRate)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarInfiltration)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarExfiltration)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarSurfaceRunoff)             = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%mLayerSatHydCondMP)              = flux2state(state1=integerMissing,  state2=integerMissing)
 flux2state_orig(iLookFLUX%mLayerSatHydCond)                = flux2state(state1=integerMissing,  state2=integerMissing)
 flux2state_orig(iLookFLUX%iLayerSatHydCond)                = flux2state(state1=integerMissing,  state2=integerMissing)
 flux2state_orig(iLookFLUX%mLayerHydCond)                   = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%iLayerLiqFluxSoil)               = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%mLayerLiqFluxSoil)               = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%mLayerBaseflow)                  = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%mLayerColumnInflow)              = flux2state(state1=integerMissing,  state2=integerMissing)
 flux2state_orig(iLookFLUX%mLayerColumnOutflow)             = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarSoilBaseflow)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarSoilDrainage)              = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarAquiferRecharge)           = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarAquiferTranspire)          = flux2state(state1=iname_matLayer,  state2=integerMissing)
 flux2state_orig(iLookFLUX%scalarAquiferBaseflow)           = flux2state(state1=iname_matLayer,  state2=integerMissing)

 ! ** copy across flux metadata
 do iVar=1,nFlux
  flux2state_orig(iVar)%var_info = flux_meta(iVar)
 end do

 ! ** check all variables are defined
 do iVar=1,nFlux
  if(flux2state_orig(iVar)%state1==integerUndefined .or. flux2state_orig(iVar)%state2==integerUndefined)then
   message=trim(message)//'flux-to-state mapping is undefined for variable "'//trim(flux_meta(iVar)%varname)//'"'
   err=20; return
  endif
 end do

 ! -----
 ! - liquid water state variables...
 ! ---------------------------------

 ! initialize to the original structure
 do iVar=1,nFlux
  flux2state_liq(iVar)%state1 = flux2state_orig(iVar)%state1 
  flux2state_liq(iVar)%state2 = flux2state_orig(iVar)%state2
 end do

 ! modify the state type names associated with the flux mapping structure
 do iVar=1,nFlux
  ! (mass of total water on the vegetation canopy --> mass of liquid water)
  if(flux2state_liq(iVar)%state1==iname_watCanopy) flux2state_liq(iVar)%state1=iname_liqCanopy
  if(flux2state_liq(iVar)%state2==iname_watCanopy) flux2state_liq(iVar)%state2=iname_liqCanopy
  ! (volumetric total water in the snow+soil domain --> volumetric liquid water)
  if(flux2state_liq(iVar)%state1==iname_watLayer)  flux2state_liq(iVar)%state1=iname_liqLayer
  if(flux2state_liq(iVar)%state2==iname_watLayer)  flux2state_liq(iVar)%state2=iname_liqLayer
  ! (total water matric potential in the snow+soil domain --> liquid water matric potential)
  if(flux2state_liq(iVar)%state1==iname_matLayer)  flux2state_liq(iVar)%state1=iname_lmpLayer
  if(flux2state_liq(iVar)%state2==iname_matLayer)  flux2state_liq(iVar)%state2=iname_lmpLayer
 end do

 ! copy across flux metadata
 do iVar=1,nFlux
  flux2state_liq(iVar)%var_info = flux_meta(iVar)
 end do

 end subroutine flxMapping

end module flxMapping_module
