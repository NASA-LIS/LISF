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

module diagn_evar_module
! data types
USE nrtype

! physical constants
USE multiconst,only:&
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  & ! intrinsic density of water    (kg m-3)
                    ! specific heat
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_ice,      & ! specific heat of ice          (J kg-1 K-1)
                    Cp_soil,     & ! specific heat of soil         (J kg-1 K-1)
                    Cp_water,    & ! specific heat of liquid water (J kg-1 K-1)
                    ! thermal conductivity
                    lambda_air,  & ! thermal conductivity of air   (J s-1 m-1)
                    lambda_ice,  & ! thermal conductivity of ice   (J s-1 m-1)
                    lambda_water   ! thermal conductivity of water (J s-1 m-1)

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables that define the layer type
USE globalData,only:iname_snow     ! snow
USE globalData,only:iname_soil     ! soil

! model decisions
USE mDecisions_module,only:Smirnova2000  ! option for temporally constant thermal conductivity

implicit none
private
public::diagn_evar
! algorithmic parameters
real(dp),parameter     :: valueMissing=-9999._dp  ! missing value, used when diagnostic or state variables are undefined
real(dp),parameter     :: verySmall=1.e-6_dp   ! used as an additive constant to check if substantial difference among real numbers
real(dp),parameter     :: mpe=1.e-6_dp         ! prevents overflow error if division by zero
real(dp),parameter     :: dx=1.e-6_dp          ! finite difference increment
contains


 ! **********************************************************************************************************
 ! public subroutine diagn_evar: compute diagnostic energy variables (thermal conductivity and heat capacity)
 ! **********************************************************************************************************
 subroutine diagn_evar(&
                       ! input: control variables
                       computeVegFlux,          & ! intent(in): flag to denote if computing the vegetation flux
                       canopyDepth,             & ! intent(in): canopy depth (m)
                       ! input/output: data structures
                       mpar_data,               & ! intent(in):    model parameters
                       indx_data,               & ! intent(in):    model layer indices
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       ! output: error control
                       err,message)               ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_types,only:&
                     var_d,            & ! data vector (dp)
                     var_ilength,      & ! data vector with variable length dimension (i4b)
                     var_dlength         ! data vector with variable length dimension (dp)
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookPARAM,iLookPROG,iLookDIAG,iLookINDEX  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
 ! provide access to named variables for thermal conductivity of soil
 USE globalData,only:model_decisions        ! model decision structure
 USE mDecisions_module,only: funcSoilWet, & ! function of soil wetness
                             mixConstit,  & ! mixture of constituents
                             hanssonVZJ     ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
 ! provide access to external subroutines
 USE snow_utils_module,only:tcond_snow            ! compute thermal conductivity of snow
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)         :: computeVegFlux         ! logical flag to denote if computing the vegetation flux
 real(dp),intent(in)             :: canopyDepth            ! depth of the vegetation canopy (m)
 ! input/output: data structures
 type(var_dlength),intent(in)    :: mpar_data              ! model parameters
 type(var_ilength),intent(in)    :: indx_data              ! model layer indices
 type(var_dlength),intent(in)    :: prog_data              ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data              ! model diagnostic variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)                :: cmessage               ! error message of downwind routine
 integer(i4b)                      :: iLayer                 ! index of model layer
 integer(i4b)                      :: iSoil                  ! index of soil layer
 real(dp)                          :: TCn                    ! thermal conductivity below the layer interface (W m-1 K-1)
 real(dp)                          :: TCp                    ! thermal conductivity above the layer interface (W m-1 K-1)
 real(dp)                          :: zdn                    ! height difference between interface and lower value (m)
 real(dp)                          :: zdp                    ! height difference between interface and upper value (m)
 real(dp)                          :: bulkden_soil           ! bulk density of soil (kg m-3)
 real(dp)                          :: lambda_drysoil         ! thermal conductivity of dry soil (W m-1)
 real(dp)                          :: lambda_wetsoil         ! thermal conductivity of wet soil (W m-1)
 real(dp)                          :: lambda_wet             ! thermal conductivity of the wet material
 real(dp)                          :: relativeSat            ! relative saturation (-)
 real(dp)                          :: kerstenNum             ! the Kersten number (-), defining weight applied to conductivity of the wet medium
 real(dp)                          :: den                    ! denominator in the thermal conductivity calculations
 ! local variables to reproduce the thermal conductivity of Hansson et al. VZJ 2005
 real(dp),parameter                :: c1=0.55_dp             ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(dp),parameter                :: c2=0.8_dp              ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(dp),parameter                :: c3=3.07_dp             ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(dp),parameter                :: c4=0.13_dp             ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(dp),parameter                :: c5=4._dp               ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(dp),parameter                :: f1=13.05_dp            ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(dp),parameter                :: f2=1.06_dp             ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(dp)                          :: fArg,xArg              ! temporary variables (see Hansson et al. VZJ 2005 for details)
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in data structure
 associate(&
 ! input: model decisions
 ixThCondSnow            => model_decisions(iLookDECISIONS%thCondSnow)%iDecision,      & ! intent(in): choice of method for thermal conductivity of snow
 ixThCondSoil            => model_decisions(iLookDECISIONS%thCondSoil)%iDecision,      & ! intent(in): choice of method for thermal conductivity of soil
 ! input: state variables
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),           & ! intent(in): canopy ice content (kg m-2)
 scalarCanopyLiquid      => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),           & ! intent(in): canopy liquid water content (kg m-2)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat,             & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat,             & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
 ! input: coordinate variables
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1),                    & ! intent(in): number of snow layers 
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1),                    & ! intent(in): number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1),                  & ! intent(in): total number of layers
 layerType               => indx_data%var(iLookINDEX%layerType)%dat,                   & ! intent(in): layer type (iname_soil or iname_snow)
 mLayerHeight            => prog_data%var(iLookPROG%mLayerHeight)%dat,                 & ! intent(in): height at the mid-point of each layer (m)
 iLayerHeight            => prog_data%var(iLookPROG%iLayerHeight)%dat,                 & ! intent(in): height at the interface of each layer (m)
 ! input: heat capacity and thermal conductivity
 specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),          & ! intent(in): specific heat of vegetation (J kg-1 K-1)
 maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),        & ! intent(in): maximum mass of vegetation (kg m-2)
 fixedThermalCond_snow   => mpar_data%var(iLookPARAM%fixedThermalCond_snow)%dat(1),    & ! intent(in): temporally constant thermal conductivity of snow (W m-1 K-1)
 ! input: depth varying soil parameters
 iden_soil               => mpar_data%var(iLookPARAM%soil_dens_intr)%dat,              & ! intent(in): intrinsic density of soil (kg m-3)
 thCond_soil             => mpar_data%var(iLookPARAM%thCond_soil)%dat,                 & ! intent(in): thermal conductivity of soil (W m-1 K-1)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat,                   & ! intent(in): soil porosity (-)
 frac_sand               => mpar_data%var(iLookPARAM%frac_sand)%dat,                   & ! intent(in): fraction of sand (-)
 frac_silt               => mpar_data%var(iLookPARAM%frac_silt)%dat,                   & ! intent(in): fraction of silt (-)
 frac_clay               => mpar_data%var(iLookPARAM%frac_clay)%dat,                   & ! intent(in): fraction of clay (-)
 ! output: diagnostic variables
 scalarBulkVolHeatCapVeg => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),   & ! intent(out): volumetric heat capacity of the vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat,           & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1)
 mLayerThermalC          => diag_data%var(iLookDIAG%mLayerThermalC)%dat,               & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
 iLayerThermalC          => diag_data%var(iLookDIAG%iLayerThermalC)%dat,               & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolFracAir        => diag_data%var(iLookDIAG%mLayerVolFracAir)%dat              & ! intent(out): volumetric fraction of air in each layer (-)
 )  ! end associate statement
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="diagn_evar/"

 ! initialize the soil layer
 iSoil=integerMissing

 ! compute the bulk volumetric heat capacity of vegetation (J m-3 K-1)
 if(computeVegFlux)then
  scalarBulkVolHeatCapVeg = specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                            Cp_water*scalarCanopyLiquid/canopyDepth       + & ! liquid water component
                            Cp_ice*scalarCanopyIce/canopyDepth                ! ice component
 else
  scalarBulkVolHeatCapVeg = valueMissing
 end if
 !print*, 'diagn_evar: scalarBulkVolHeatCapVeg = ', scalarBulkVolHeatCapVeg

 ! loop through layers
 do iLayer=1,nLayers

  ! get the soil layer
  if(iLayer>nSnow) iSoil = iLayer-nSnow

  ! compute the thermal conductivity of dry and wet soils (W m-1)
  ! NOTE: this is actually constant over the simulation, and included here for clarity
  if(ixThCondSoil == funcSoilWet .and. layerType(iLayer)==iname_soil)then
   bulkden_soil   = iden_soil(iSoil)*( 1._dp - theta_sat(iSoil) )
   lambda_drysoil = (0.135_dp*bulkden_soil + 64.7_dp) / (iden_soil(iSoil) - 0.947_dp*bulkden_soil)
   lambda_wetsoil = (8.80_dp*frac_sand(iSoil) + 2.92_dp*frac_clay(iSoil)) / (frac_sand(iSoil) + frac_clay(iSoil))
  end if

  ! *****
  ! * compute the volumetric fraction of air in each layer...
  ! *********************************************************
  select case(layerType(iLayer))
   case(iname_soil); mLayerVolFracAir(iLayer) = theta_sat(iSoil) - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))
   case(iname_snow); mLayerVolFracAir(iLayer) = 1._dp - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))
   case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute volumetric fraction of air'; return
  end select

  ! *****
  ! * compute the volumetric heat capacity of each layer (J m-3 K-1)...
  ! *******************************************************************
  select case(layerType(iLayer))
   ! * soil
   case(iname_soil)
    mLayerVolHtCapBulk(iLayer) = iden_soil(iSoil)  * Cp_soil  * ( 1._dp - theta_sat(iSoil) ) + & ! soil component
                                 iden_ice          * Cp_Ice   * mLayerVolFracIce(iLayer)     + & ! ice component
                                 iden_water        * Cp_water * mLayerVolFracLiq(iLayer)     + & ! liquid water component
                                 iden_air          * Cp_air   * mLayerVolFracAir(iLayer)         ! air component
   ! * snow
   case(iname_snow)
    mLayerVolHtCapBulk(iLayer) = iden_ice          * Cp_ice   * mLayerVolFracIce(iLayer)     + & ! ice component
                                 iden_water        * Cp_water * mLayerVolFracLiq(iLayer)     + & ! liquid water component
                                 iden_air          * Cp_air   * mLayerVolFracAir(iLayer)         ! air component
   case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute olumetric heat capacity'; return
  end select

  ! *****
  ! * compute the thermal conductivity of snow and soil at the mid-point of each layer...
  ! *************************************************************************************
  select case(layerType(iLayer))

   ! ***** soil
   case(iname_soil)

    ! select option for thermal conductivity of soil
    select case(ixThCondSoil)

     ! ** function of soil wetness
     case(funcSoilWet)

      ! compute the thermal conductivity of the wet material (W m-1)
      lambda_wet  = lambda_wetsoil**( 1._dp - theta_sat(iSoil) ) * lambda_water**theta_sat(iSoil) * lambda_ice**(theta_sat(iSoil) - mLayerVolFracLiq(iLayer))
      relativeSat = (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))/theta_sat(iSoil)  ! relative saturation
      ! compute the Kersten number (-)
      if(relativeSat > 0.1_dp)then ! log10(0.1) = -1
       kerstenNum = log10(relativeSat) + 1._dp
      else
       kerstenNum = 0._dp  ! dry thermal conductivity
      endif
      ! ...and, compute the thermal conductivity
      mLayerThermalC(iLayer) = kerstenNum*lambda_wet + (1._dp - kerstenNum)*lambda_drysoil

     ! ** mixture of constituents
     case(mixConstit)
      mLayerThermalC(iLayer) = thCond_soil(iSoil) * ( 1._dp - theta_sat(iSoil) ) + & ! soil component
                               lambda_ice         * mLayerVolFracIce(iLayer)     + & ! ice component
                               lambda_water       * mLayerVolFracLiq(iLayer)     + & ! liquid water component
                               lambda_air         * mLayerVolFracAir(iLayer)         ! air component

     ! ** test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
     case(hanssonVZJ)
      fArg  = 1._dp + f1*mLayerVolFracIce(iLayer)**f2
      xArg  = mLayerVolFracLiq(iLayer) + fArg*mLayerVolFracIce(iLayer)
      mLayerThermalC(iLayer) = c1 + c2*xArg + (c1 - c4)*exp(-(c3*xArg)**c5)

     ! ** check
     case default; err=20; message=trim(message)//'unable to identify option for thermal conductivity of soil'; return

    end select  ! option for the thermal conductivity of soil
    
   ! ***** snow
   case(iname_snow)
    ! temporally constant thermal conductivity
    if(ixThCondSnow==Smirnova2000)then
     mLayerThermalC(iLayer) = fixedThermalCond_snow
    ! thermal conductivity as a function of snow density
    else
     call tcond_snow(mLayerVolFracIce(iLayer)*iden_ice,  & ! input: snow density (kg m-3)
                     mLayerThermalC(iLayer),             & ! output: thermal conductivity (W m-1 K-1)
                     err,cmessage)                         ! output: error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    endif

   ! * error check
   case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute thermal conductivity'; return

  end select
  !print*, 'iLayer, mLayerThermalC(iLayer) = ', iLayer, mLayerThermalC(iLayer)

 end do  ! looping through layers
 !pause

 ! *****
 ! * compute the thermal conductivity of snow at the interface of each layer...
 ! ****************************************************************************
 do iLayer=1,nLayers-1  ! (loop through layers)
  ! get temporary variables
  TCn = mLayerThermalC(iLayer)    ! thermal conductivity below the layer interface (W m-1 K-1)
  TCp = mLayerThermalC(iLayer+1)  ! thermal conductivity above the layer interface (W m-1 K-1)
  zdn = iLayerHeight(iLayer)   - mLayerHeight(iLayer) ! height difference between interface and lower value (m)
  zdp = mLayerHeight(iLayer+1) - iLayerHeight(iLayer) ! height difference between interface and upper value (m)
  den = TCn*zdp + TCp*zdn  ! denominator
  ! compute thermal conductivity
  if(TCn+TCp > epsilon(TCn))then
   iLayerThermalC(iLayer) = (TCn*TCp*(zdn + zdp)) / den
  else
   iLayerThermalC(iLayer) = (TCn*zdn +  TCp*zdp) / (zdn + zdp)
  endif
  !write(*,'(a,1x,i4,1x,10(f9.3,1x))') 'iLayer, TCn, TCp, zdn, zdp, iLayerThermalC(iLayer) = ', iLayer, TCn, TCp, zdn, zdp, iLayerThermalC(iLayer)
 end do  ! looping through layers

 ! special case of hansson
 if(ixThCondSoil==hanssonVZJ)then
  iLayerThermalC(0) = 28._dp*(0.5_dp*(iLayerHeight(1) - iLayerHeight(0)))
 else
  iLayerThermalC(0) = mLayerThermalC(1)
 end if

 ! assume the thermal conductivity at the domain boundaries is equal to the thermal conductivity of the layer
 iLayerThermalC(nLayers) = mLayerThermalC(nLayers)

 ! end association to variables in the data structure
 end associate

 end subroutine diagn_evar


end module diagn_evar_module
