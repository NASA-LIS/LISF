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

module check_icond_module
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number

! define modeling decisions
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation

implicit none
private
public::check_icond
contains

 ! ************************************************************************************************
 ! public subroutine check_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine check_icond(nGRU,                          & ! number of GRUs and HRUs
                        progData,                      & ! model prognostic (state) variables
                        mparData,                      & ! model parameters
                        indxData,                      & ! layer index data
                        err,message)                     ! error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookParam                          ! variable lookup structure
 USE var_lookup,only:iLookProg                           ! variable lookup structure
 USE var_lookup,only:iLookIndex                          ! variable lookup structure
 USE globalData,only:gru_struc                           ! gru-hru mapping structures
 USE data_types,only:gru_hru_doubleVec                   ! actual data
 USE data_types,only:gru_hru_intVec                      ! actual data
 USE globaldata,only:iname_soil,iname_snow               ! named variables to describe the type of layer
 USE multiconst,only:&
                       LH_fus,    &                      ! latent heat of fusion                (J kg-1)
                       iden_ice,  &                      ! intrinsic density of ice             (kg m-3)
                       iden_water,&                      ! intrinsic density of liquid water    (kg m-3)
                       gravity,   &                      ! gravitational acceleration           (m s-2)
                       Tfreeze                           ! freezing point of pure water         (K)
 USE snow_utils_module,only:fracliquid                   ! compute volumetric fraction of liquid water in snow based on temperature
 USE updatState_module,only:updateSnow                   ! update snow states
 USE updatState_module,only:updateSoil                   ! update soil states
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 integer(i4b)           ,intent(in)    :: nGRU           ! number of grouped response units
 type(gru_hru_doubleVec),intent(inout) :: progData       ! prognostic vars 
 type(gru_hru_doubleVec),intent(in)    :: mparData       ! parameters 
 type(gru_hru_intVec)   ,intent(in)    :: indxData       ! layer indexes 
 integer(i4b)           ,intent(out)   :: err            ! error code
 character(*)           ,intent(out)   :: message        ! returned error message

 ! locals
 character(len=256)                     :: cmessage      ! downstream error message
 integer(i4b)                           :: iGRU          ! loop index 
 integer(i4b)                           :: iHRU          ! loop index 

 ! temporary variables for realism checks
 integer(i4b)                      :: iLayer             ! index of model layer
 integer(i4b)                      :: iSoil              ! index of soil layer
 real(dp)                          :: fLiq               ! fraction of liquid water on the vegetation canopy (-)
 real(dp)                          :: vGn_m              ! van Genutchen "m" parameter (-)
 real(dp)                          :: tWat               ! total water on the vegetation canopy (kg m-2)
 real(dp)                          :: scalarTheta        ! liquid water equivalent of total water [liquid water + ice] (-)
 real(dp)                          :: h1,h2              ! used to check depth and height are consistent
 integer(i4b)                      :: nLayers            ! total number of layers
 real(dp)                          :: kappa              ! constant in the freezing curve function (m K-1)
 integer(i4b)                      :: nSnow              ! number of snow layers

 ! --------------------------------------------------------------------------------------------------------

 ! Start procedure here
 err=0; message="check_icond/"

 ! --------------------------------------------------------------------------------------------------------
 ! Check that the initial conditions do not conflict with parameters, structure, etc.
 ! --------------------------------------------------------------------------------------------------------
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   ! ensure the spectral average albedo is realistic
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) > mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMax)%dat(1)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMax)%dat(1)
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) < mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinWinter)%dat(1)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinWinter)%dat(1)
   ! ensure the visible albedo is realistic
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) > mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMaxVisible)%dat(1)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMaxVisible)%dat(1)
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) < mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinVisible)%dat(1)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinVisible)%dat(1)
   ! ensure the nearIR albedo is realistic
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) > mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMaxNearIR)%dat(1)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMaxNearIR)%dat(1)
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) < mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinNearIR)%dat(1)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinNearIR)%dat(1)
  end do
 end do
 
 ! ensure the initial conditions are consistent with the constitutive functions
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
 
   ! associate local variables with variables in the data structures
   associate(&
   ! state variables in the vegetation canopy
   scalarCanopyTemp  => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyTemp)%dat(1)   , & ! canopy temperature
   scalarCanopyIce   => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyIce)%dat(1)    , & ! mass of ice on the vegetation canopy (kg m-2)
   scalarCanopyLiq   => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyLiq)%dat(1)    , & ! mass of liquid water on the vegetation canopy (kg m-2)
   ! state variables in the snow+soil domain
   mLayerTemp        => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerTemp)%dat            , & ! temperature (K)
   mLayerVolFracLiq  => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracLiq)%dat      , & ! volumetric fraction of liquid water in each snow layer (-)
   mLayerVolFracIce  => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracIce)%dat      , & ! volumetric fraction of ice in each snow layer (-)
   mLayerMatricHead  => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerMatricHead)%dat      , & ! matric head (m)
   mLayerLayerType   => indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat            , & ! type of layer (ix_soil or ix_snow)
   ! depth varying soil properties
   vGn_alpha         => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%vGn_alpha)%dat            , & ! van Genutchen "alpha" parameter (m-1)
   vGn_n             => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%vGn_n)%dat                , & ! van Genutchen "n" parameter (-)
   theta_sat         => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%theta_sat)%dat            , & ! soil porosity (-)
   theta_res         => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%theta_res)%dat            , & ! soil residual volumetric water content (-)
   ! snow parameters
   snowfrz_scale     => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%snowfrz_scale)%dat(1)     , & ! scaling parameter for the snow freezing curve (K-1)
   FCapil            => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%FCapil)%dat(1)              & ! fraction of pore space in tension storage (-)
   )  ! (associate local variables with model parameters)

   ! compute the constant in the freezing curve function (m K-1)
   kappa  = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2

   ! modify the liquid water and ice in the canopy
   if(scalarCanopyIce > 0._dp .and. scalarCanopyTemp > Tfreeze)then
    message=trim(message)//'canopy ice > 0 when canopy temperature > Tfreeze'
    err=20; return
   end if
   fLiq = fracliquid(scalarCanopyTemp,snowfrz_scale)  ! fraction of liquid water (-)
   tWat = scalarCanopyLiq + scalarCanopyIce           ! total water (kg m-2)
   scalarCanopyLiq = fLiq*tWat                        ! mass of liquid water on the canopy (kg m-2)
   scalarCanopyIce = (1._dp - fLiq)*tWat              ! mass of ice on the canopy (kg m-2)

   ! number of layers
   nLayers = gru_struc(iGRU)%hruInfo(iHRU)%nSnow + gru_struc(iGRU)%hruInfo(iHRU)%nSoil
   nSnow   = gru_struc(iGRU)%hruInfo(iHRU)%nSnow

   ! loop through all layers
   do iLayer=1,nLayers

    ! compute liquid water equivalent of total water (liquid plus ice)
    if (iLayer>nSnow) then ! soil layer = no volume expansion
     iSoil       = iLayer - nSnow
     vGn_m       = 1._dp - 1._dp/vGn_n(iSoil)
     scalarTheta = mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)
    else ! snow layer = volume expansion allowed
     iSoil       = integerMissing
     vGn_m       = realMissing
     scalarTheta = mLayerVolFracIce(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiq(iLayer)
    end if

    ! *****
    ! * check that the initial volumetric fraction of liquid water and ice is reasonable...
    ! *************************************************************************************
    select case(mlayerLayerType(iLayer))

     ! ***** snow
     case(iname_snow)
      ! (check liquid water)
      if(mLayerVolFracLiq(iLayer) < 0._dp)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water < 0: layer = ',iLayer; err=20; return; end if
      if(mLayerVolFracLiq(iLayer) > 1._dp)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water > 1: layer = ',iLayer; err=20; return; end if
      ! (check ice)
      if(mLayerVolFracIce(iLayer) > 0.80_dp)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice > 0.80: layer = ',iLayer; err=20; return; end if
      if(mLayerVolFracIce(iLayer) < 0.05_dp)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0.05: layer = ',iLayer; err=20; return; end if
      ! check total water
      if(scalarTheta > 0.80_dp)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] > 0.80: layer = ',iLayer; err=20; return; end if
      if(scalarTheta < 0.05_dp)then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] < 0.05: layer = ',iLayer; err=20; return; end if

     ! ***** soil
     case(iname_soil)
      ! (check liquid water)
      if(mLayerVolFracLiq(iLayer) < theta_res(iSoil) )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water < theta_res: layer = ',iLayer; err=20; return; end if
      if(mLayerVolFracLiq(iLayer) > theta_sat(iSoil) )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water > theta_sat: layer = ',iLayer; err=20; return; end if
      ! (check ice)
      if(mLayerVolFracIce(iLayer) < 0._dp            )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0: layer = '        ,iLayer; err=20; return; end if
      if(mLayerVolFracIce(iLayer) > theta_sat(iSoil) )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice > theta_sat: layer = ',iLayer; err=20; return; end if
      ! check total water
      if(scalarTheta < theta_res(iSoil) )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] < theta_res: layer = ',iLayer; err=20; return; end if 
      if(scalarTheta > theta_sat(iSoil) )then; write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with total water fraction [liquid + ice] > theta_sat: layer = ',iLayer; err=20; return; end if

     case default
      err=20; message=trim(message)//'cannot identify layer type'; return
    end select

    ! *****
    ! * check that the initial conditions are consistent with the constitutive functions...
    ! *************************************************************************************
    select case(mLayerLayerType(iLayer))
  
     ! ** snow
     case(iname_snow)
  
      ! check that snow temperature is less than freezing
      if(mLayerTemp(iLayer) > Tfreeze)then
       message=trim(message)//'initial snow temperature is greater than freezing'
       err=20; return
      end if
  
      ! ensure consistency among state variables
      call updateSnow(&
                      ! input
                      mLayerTemp(iLayer),                        & ! intent(in): temperature (K)
                      scalarTheta,                               & ! intent(in): mass fraction of total water (-)
                      snowfrz_scale,                             & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                      ! output
                      mLayerVolFracLiq(iLayer),                  & ! intent(out): volumetric fraction of liquid water (-)
                      mLayerVolFracIce(iLayer),                  & ! intent(out): volumetric fraction of ice (-)
                      fLiq,                                      & ! intent(out): fraction of liquid water (-)
                      err,cmessage)                                ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
  
     ! ** soil
     case(iname_soil)
  
      ! ensure consistency among state variables
      call updateSoil(&
                      ! input
                      mLayerTemp(iLayer),                                                    & ! intent(in): layer temperature (K)
                      mLayerMatricHead(iLayer-nSnow),                                        & ! intent(in): matric head (m)
                      vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m, & ! intent(in): van Genutchen soil parameters
                      ! output
                      scalarTheta,                                                           & ! intent(out): volumetric fraction of total water (-)
                      mLayerVolFracLiq(iLayer),                                              & ! intent(out): volumetric fraction of liquid water (-)
                      mLayerVolFracIce(iLayer),                                              & ! intent(out): volumetric fraction of ice (-)
                      err,cmessage)                                                            ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
  
     case default; err=10; message=trim(message)//'unknown case for model layer'; return
    end select
  
   end do  ! (looping through layers)
  
   ! end association to variables in the data structures
   end associate
  
   ! if snow layers exist, compute snow depth and SWE
   if(nSnow > 0)then
    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSWE)%dat(1)       = sum( (progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                                               progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice)        * &
                                                                               progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
   end if  ! if snow layers exist

   ! check that the layering is consistent
   do iLayer=1,nLayers
    h1 = sum(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(1:iLayer)) ! sum of the depths up to the current layer
    h2 = progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%iLayerHeight)%dat(iLayer) - progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%iLayerHeight)%dat(0)  ! difference between snow-atm interface and bottom of layer
    if(abs(h1 - h2) > 1.e-12_dp)then
     write(message,'(a,1x,i0)') trim(message)//'mis-match between layer depth and layer height [suggest round numbers in initial conditions file]; layer = ', iLayer
     err=20; return
    end if
   end do

  end do ! iHRU
 end do ! iGRU

 end subroutine check_icond

end module check_icond_module
