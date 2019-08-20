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

module volicePack_module
! numerical recipes data types
USE nrtype
! physical constants
USE multiconst,only:&
                    Tfreeze,  & ! freezing point              (K)
                    LH_fus,   & ! latent heat of fusion       (J kg-1)
                    LH_vap,   & ! latent heat of vaporization (J kg-1)
                    LH_sub,   & ! latent heat of sublimation  (J kg-1)
                    iden_air, & ! intrinsic density of air    (kg m-3)
                    iden_ice, & ! intrinsic density of ice    (kg m-3)
                    iden_water  ! intrinsic density of water  (kg m-3)
implicit none
private
public::volicePack
public::newsnwfall

contains


 ! ************************************************************************************************
 ! public subroutine volicePack: combine and sub-divide layers if necessary)
 ! ************************************************************************************************
 subroutine volicePack(&
                       ! input/output: model data structures
                       tooMuchMelt,                 & ! intent(in):    flag to force merge of snow layers
                       model_decisions,             & ! intent(in):    model decisions
                       mpar_data,                   & ! intent(in):    model parameters
                       indx_data,                   & ! intent(inout): type of each layer
                       prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,                   & ! intent(inout): model fluxes for a local HRU
                       ! output
                       modifiedLayers,              & ! intent(out): flag to denote that layers were modified
                       err,message)                   ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_types,only:&
                     var_d,            & ! data vector (dp)
                     var_ilength,      & ! data vector with variable length dimension (i4b)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! external subroutine
 USE layerMerge_module,only:layerMerge   ! merge snow layers if they are too thin
 USE layerDivide_module,only:layerDivide ! sub-divide layers if they are too thick
 implicit none
 ! ------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 logical(lgt),intent(in)         :: tooMuchMelt         ! flag to denote that ice is insufficient to support melt
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 type(var_dlength),intent(in)    :: mpar_data           ! model parameters
 type(var_ilength),intent(inout) :: indx_data           ! type of each layer
 type(var_dlength),intent(inout) :: prog_data           ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data           ! model diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data           ! model flux variables
 ! output
 logical(lgt),intent(out)        :: modifiedLayers      ! flag to denote that we modified the layers
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! ------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)              :: cmessage            ! error message of downwind routine
 logical(lgt)                    :: mergedLayers        ! flag to denote that layers were merged
 logical(lgt)                    :: divideLayer         ! flag to denote that a layer was divided
 ! initialize error control
 err=0; message='volicePack/'

 ! divide snow layers if too thick
 call layerDivide(&
                  ! input/output: model data structures
                  model_decisions,             & ! intent(in):    model decisions
                  mpar_data,                   & ! intent(in):    model parameters
                  indx_data,                   & ! intent(inout): type of each layer
                  prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                  diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                  flux_data,                   & ! intent(inout): model fluxes for a local HRU
                  ! output
                  divideLayer,                 & ! intent(out): flag to denote that layers were modified
                  err,cmessage)                  ! intent(out): error control
 if(err/=0)then; err=65; message=trim(message)//trim(cmessage); return; end if

 ! merge snow layers if they are too thin
 call layerMerge(&
                 ! input/output: model data structures
                 tooMuchMelt,                 & ! intent(in):    flag to force merge of snow layers
                 model_decisions,             & ! intent(in):    model decisions
                 mpar_data,                   & ! intent(in):    model parameters
                 indx_data,                   & ! intent(inout): type of each layer
                 prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                 diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                 flux_data,                   & ! intent(inout): model fluxes for a local HRU
                 ! output
                 mergedLayers,                & ! intent(out): flag to denote that layers were modified
                 err,cmessage)                  ! intent(out): error control
 if(err/=0)then; err=65; message=trim(message)//trim(cmessage); return; end if

 ! flag if layers were modified
 modifiedLayers = (mergedLayers .or. divideLayer)

 end subroutine volicePack


 ! ************************************************************************************************
 ! public subroutine newsnwfall: add new snowfall to the system
 ! ************************************************************************************************
 subroutine newsnwfall(&
                       ! input: model control
                       dt,                        & ! time step (seconds)
                       snowLayers,                & ! logical flag if snow layers exist
                       fc_param,                  & ! freeezing curve parameter for snow (K-1)
                       ! input: diagnostic scalar variables
                       scalarSnowfallTemp,        & ! computed temperature of fresh snow (K)
                       scalarNewSnowDensity,      & ! computed density of new snow (kg m-3)
                       scalarThroughfallSnow,     & ! throughfall of snow through the canopy (kg m-2 s-1)
                       scalarCanopySnowUnloading, & ! unloading of snow from the canopy (kg m-2 s-1)
                       ! input/output: state variables
                       scalarSWE,                 & ! SWE (kg m-2)
                       scalarSnowDepth,           & ! total snow depth (m)
                       surfaceLayerTemp,          & ! temperature of each layer (K)
                       surfaceLayerDepth,         & ! depth of each layer (m)
                       surfaceLayerVolFracIce,    & ! volumetric fraction of ice in each layer (-)
                       surfaceLayerVolFracLiq,    & ! volumetric fraction of liquid water in each layer (-)
                       ! output: error control
                       err,message                ) ! error control
 ! computational modules
 USE snow_utils_module,only:fracliquid,templiquid                  ! functions to compute temperature/liquid water
 ! add new snowfall to the system
 implicit none
 ! input: model control
 real(dp),intent(in)                 :: dt                         ! time step (seconds)
 logical(lgt),intent(in)             :: snowLayers                 ! logical flag if snow layers exist
 real(dp),intent(in)                 :: fc_param                   ! freeezing curve parameter for snow (K-1)
 ! input: diagnostic scalar variables
 real(dp),intent(in)                 :: scalarSnowfallTemp         ! computed temperature of fresh snow (K)
 real(dp),intent(in)                 :: scalarNewSnowDensity       ! computed density of new snow (kg m-3)
 real(dp),intent(in)                 :: scalarThroughfallSnow      ! throughfall of snow through the canopy (kg m-2 s-1)
 real(dp),intent(in)                 :: scalarCanopySnowUnloading  ! unloading of snow from the canopy (kg m-2 s-1)
 ! input/output: state variables
 real(dp),intent(inout)              :: scalarSWE                  ! SWE (kg m-2)
 real(dp),intent(inout)              :: scalarSnowDepth            ! total snow depth (m)
 real(dp),intent(inout)              :: surfaceLayerTemp           ! temperature of each layer (K)
 real(dp),intent(inout)              :: surfaceLayerDepth          ! depth of each layer (m)
 real(dp),intent(inout)              :: surfaceLayerVolFracIce     ! volumetric fraction of ice in each layer (-)
 real(dp),intent(inout)              :: surfaceLayerVolFracLiq     ! volumetric fraction of liquid water in each layer (-)
 ! output: error control
 integer(i4b),intent(out)            :: err                        ! error code
 character(*),intent(out)            :: message                    ! error message
 ! define local variables
 real(dp)                            :: newSnowfall                ! new snowfall -- throughfall and unloading (kg m-2 s-1)
 real(dp)                            :: newSnowDepth               ! new snow depth (m)
 real(dp),parameter                  :: densityCanopySnow=200._dp  ! density of snow on the vegetation canopy (kg m-3)
 real(dp)                            :: totalMassIceSurfLayer      ! total mass of ice in the surface layer (kg m-2)
 real(dp)                            :: totalDepthSurfLayer        ! total depth of the surface layer (m)
 real(dp)                            :: volFracWater               ! volumetric fraction of total water, liquid and ice (-)
 real(dp)                            :: fracLiq                    ! fraction of liquid water (-)
 real(dp)                            :: SWE                        ! snow water equivalent after snowfall (kg m-2)
 real(dp)                            :: tempSWE0                   ! temporary SWE before snowfall, used to check mass balance (kg m-2)
 real(dp)                            :: tempSWE1                   ! temporary SWE after snowfall, used to check mass balance (kg m-2)
 real(dp)                            :: xMassBalance               ! mass balance check (kg m-2)
 real(dp),parameter                  :: verySmall=1.e-8_dp         ! a very small number -- used to check mass balance
 ! initialize error control
 err=0; message="newsnwfall/"

 ! compute the new snowfall (kg m-2 s-1)
 newSnowfall = scalarThroughfallSnow + scalarCanopySnowUnloading

 ! early return if there is no snowfall
 if(newSnowfall < tiny(dt)) return

 ! compute depth of new snow
 newSnowDepth     = dt*(scalarThroughfallSnow/scalarNewSnowDensity + scalarCanopySnowUnloading/densityCanopySnow)  ! new snow depth (m)

 ! process special case of "snow without a layer"
 if(.not.snowLayers)then
  ! increment depth and water equivalent
  scalarSnowDepth = scalarSnowDepth + newSnowDepth
  scalarSWE       = scalarSWE + dt*newSnowfall

 ! add snow to the top layer (more typical case where snow layers already exist)
 else

  ! get SWE in the upper layer (used to check mass balance)
  tempSWE0 = (surfaceLayerVolFracIce*iden_ice + surfaceLayerVolFracLiq*iden_water)*surfaceLayerDepth

  ! get the total mass of liquid water and ice (kg m-2)
  totalMassIceSurfLayer  = iden_ice*surfaceLayerVolFracIce*surfaceLayerDepth + newSnowfall*dt
  ! get the total snow depth
  totalDepthSurfLayer    = surfaceLayerDepth + newSnowDepth
  !write(*,'(a,1x,10(f20.10,1x))') 'scalarSnowfallTemp, surfaceLayerTemp, newSnowDepth, surfaceLayerDepth, tempSWE0, totalMassIceSurfLayer/totalDepthSurfLayer = ', &
  !                                 scalarSnowfallTemp, surfaceLayerTemp, newSnowDepth, surfaceLayerDepth, tempSWE0, totalMassIceSurfLayer/totalDepthSurfLayer

  ! compute the new temperature
  surfaceLayerTemp       = (surfaceLayerTemp*surfaceLayerDepth + scalarSnowfallTemp*newSnowDepth) / totalDepthSurfLayer
  ! compute new SWE for the upper layer (kg m-2)
  SWE = totalMassIceSurfLayer + iden_water*surfaceLayerVolFracLiq*surfaceLayerDepth
  ! compute new volumetric fraction of liquid water and ice (-)
  volFracWater = (SWE/totalDepthSurfLayer)/iden_water
  fracLiq      = fracliquid(surfaceLayerTemp,fc_param)                           ! fraction of liquid water
  surfaceLayerVolFracIce = (1._dp - fracLiq)*volFracWater*(iden_water/iden_ice)  ! volumetric fraction of ice (-)
  surfaceLayerVolFracLiq =          fracLiq *volFracWater                        ! volumetric fraction of liquid water (-)
  ! update new layer depth (m)
  surfaceLayerDepth      = totalDepthSurfLayer

  ! get SWE in the upper layer (used to check mass balance)
  tempSWE1 = (surfaceLayerVolFracIce*iden_ice + surfaceLayerVolFracLiq*iden_water)*surfaceLayerDepth

  ! check SWE
  xMassBalance = tempSWE1 - (tempSWE0 + newSnowfall*dt)
  if (abs(xMassBalance) > verySmall)then
   write(*,'(a,1x,f20.10)') 'SWE mass balance = ', xMassBalance
   message=trim(message)//'mass balance problem'
   err=20; return
  end if

 end if  ! if snow layers already exist

 end subroutine newsnwfall


end module volicePack_module
