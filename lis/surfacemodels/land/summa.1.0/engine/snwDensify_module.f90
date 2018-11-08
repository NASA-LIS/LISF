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

module snwDensify_module
USE nrtype
USE multiconst,only:&
                    Tfreeze,    &     ! freezing point of pure water (K)
                    iden_ice,   &     ! intrinsic density of ice (kg m-3)
                    iden_water        ! intrinsic density of liquid water (kg m-3)
implicit none
private
public::snwDensify
contains

 ! ************************************************************************************************
 ! public subroutine snwDensify: compute change in snow density over the time step
 ! ************************************************************************************************
 subroutine snwDensify(&

                       ! intent(in): variables
                       dt,                             & ! intent(in): time step (s)
                       nSnow,                          & ! intent(in): number of snow layers
                       mLayerTemp,                     & ! intent(in): temperature of each layer (K)
                       mLayerMeltFreeze,               & ! intent(in): volumnetric melt in each layer (kg m-3)
                       scalarSnowSublimation,          & ! intent(in): sublimation from the snow surface (kg m-2 s-1)

                       ! intent(in): parameters
                       densScalGrowth,                 & ! intent(in): density scaling factor for grain growth (kg-1 m3)
                       tempScalGrowth,                 & ! intent(in): temperature scaling factor for grain growth (K-1)
                       grainGrowthRate,                & ! intent(in): rate of grain growth (s-1)
                       densScalOvrbdn,                 & ! intent(in): density scaling factor for overburden pressure (kg-1 m3)
                       tempScalOvrbdn,                 & ! intent(in): temperature scaling factor for overburden pressure (K-1)
                       baseViscosity,                      & ! intent(in): viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)

                       ! intent(inout): state variables
                       mLayerDepth,                    & ! intent(inout): depth of each layer (m)
                       mLayerVolFracLiqNew,            & ! intent(inout):  volumetric fraction of liquid water after itertations (-)
                       mLayerVolFracIceNew,            & ! intent(inout):  volumetric fraction of ice after itertations (-)

                       ! output: error control
                       err,message)                      ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! compute change in snow density over the time step
 implicit none
 ! intent(in): variables
 real(dp),intent(in)                 :: dt                       ! time step (seconds)
 integer(i4b),intent(in)             :: nSnow                    ! number of snow layers
 real(dp),intent(in)                 :: mLayerTemp(:)            ! temperature of each snow layer after iterations (K)
 real(dp),intent(in)                 :: mLayerMeltFreeze(:)      ! volumetric melt in each layer (kg m-3)
 real(dp),intent(in)                 :: scalarSnowSublimation    ! sublimation from the snow surface (kg m-2 s-1)
 ! intent(in): parameters
 real(dp),intent(in)                 :: densScalGrowth           ! density scaling factor for grain growth (kg-1 m3)
 real(dp),intent(in)                 :: tempScalGrowth           ! temperature scaling factor for grain growth (K-1)
 real(dp),intent(in)                 :: grainGrowthRate          ! rate of grain growth (s-1)
 real(dp),intent(in)                 :: densScalOvrbdn           ! density scaling factor for overburden pressure (kg-1 m3)
 real(dp),intent(in)                 :: tempScalOvrbdn           ! temperature scaling factor for overburden pressure (K-1)
 real(dp),intent(in)                 :: baseViscosity            ! viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
 ! intent(inout): state variables
 real(dp),intent(inout)              :: mLayerDepth(:)           ! depth of each layer (m)
 real(dp),intent(inout)              :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water in each snow layer after iterations (-)
 real(dp),intent(inout)              :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice in each snow layer after iterations (-)
 ! intent(out): error control
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! define local variables
 integer(i4b)                        :: iSnow                    ! index of snow layers
 real(dp)                            :: chi1,chi2,chi3,chi4,chi5 ! multipliers in the densification algorithm (-)
 real(dp)                            :: halfWeight               ! half of the weight of the current snow layer (kg m-2)
 real(dp)                            :: weightSnow               ! total weight of snow above the current snow layer (kg m-2)
 real(dp)                            :: CR_grainGrowth           ! compaction rate for grain growth (s-1)
 real(dp)                            :: CR_ovrvdnPress           ! compaction rate associated with over-burden pressure (s-1)
 real(dp)                            :: CR_metamorph             ! compaction rate for metamorphism (s-1)
 real(dp)                            :: massIceOld               ! mass of ice in the snow layer (kg m-2)
 real(dp)                            :: massLiqOld               ! mass of liquid water in the snow layer (kg m-2)
 real(dp)                            :: scalarDepthNew           ! updated layer depth (m)
 real(dp)                            :: volFracIceLoss           ! volumetric fraction of ice lost due to melt and sublimation (-)
 real(dp),parameter                  :: snwden_min=100._dp       ! minimum snow density for reducing metamorphism rate (kg m-3)
 real(dp),parameter                  :: snwDensityMax=550._dp    ! maximum snow density for collapse under melt (kg m-3)
 real(dp),parameter                  :: wetSnowThresh=0.01_dp    ! threshold to discriminate between "wet" and "dry" snow
 real(dp),parameter                  :: minLayerDensity=40._dp   ! minimum snow density allowed for any layer (kg m-3)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="snwDensify/"

 ! NOTE: still need to process the case of "snow without a layer"
 if(nSnow==0)return

 ! initialize the weight of snow above each layer (kg m-2)
 weightSnow = 0._dp

 ! loop through snow layers
 do iSnow=1,nSnow

  ! print starting density
  !write(*,'(a,1x,i4,1x,f9.3)') 'b4 compact: iSnow, density = ', iSnow, mLayerVolFracIceNew(iSnow)*iden_ice

  ! save mass of liquid water and ice (mass does not change)
  massIceOld = iden_ice*mLayerVolFracIceNew(iSnow)*mLayerDepth(iSnow)   ! (kg m-2)
  massLiqOld = iden_water*mLayerVolFracLiqNew(iSnow)*mLayerDepth(iSnow) ! (kg m-2)

  ! *** compute the compaction associated with grain growth (s-1)
  ! compute the base rate of grain growth (-)
  if(mLayerVolFracIceNew(iSnow)*iden_ice <snwden_min) chi1=1._dp
  if(mLayerVolFracIceNew(iSnow)*iden_ice>=snwden_min) chi1=exp(-densScalGrowth*(mLayerVolFracIceNew(iSnow)*iden_ice - snwden_min))
  ! compute the reduction of grain growth under colder snow temperatures (-)
  chi2 = exp(-tempScalGrowth*(Tfreeze - mLayerTemp(iSnow)))
  ! compute the acceleration of grain growth in the presence of liquid water (-)
  if(mLayerVolFracLiqNew(iSnow) > wetSnowThresh)then; chi3=2._dp  ! snow is "wet"
  else; chi3=1._dp; end if                                         ! snow is "dry"
  ! compute the compaction associated with grain growth (s-1)
  CR_grainGrowth = grainGrowthRate*chi1*chi2*chi3

  ! **** compute the compaction associated with over-burden pressure (s-1)
  ! compute the weight imposed on the current layer (kg m-2)
  halfWeight = (massIceOld + massLiqOld)/2._dp  ! there is some over-burden pressure from the layer itself
  weightSnow = weightSnow + halfweight          ! add half of the weight from the current layer
  ! compute the increase in compaction under colder snow temperatures (-)
  chi4 = exp(-tempScalOvrbdn*(Tfreeze - mLayerTemp(iSnow)))
  ! compute the increase in compaction under low density snow (-)
  chi5 = exp(-densScalOvrbdn*mLayerVolFracIceNew(iSnow)*iden_ice)
  ! compute the compaction associated with over-burden pressure (s-1)
  CR_ovrvdnPress = (weightSnow/baseViscosity)*chi4*chi5
  ! update the snow weight with the halfWeight not yet used
  weightSnow = weightSnow + halfweight          ! add half of the weight from the current layer

  ! *** compute the compaction rate associated with snow melt (s-1)
  ! NOTE: loss of ice due to snowmelt is implicit, so can be updated directly
  if(iden_ice*mLayerVolFracIceNew(iSnow) < snwDensityMax)then ! only collapse layers if below a critical density
   ! (compute volumetric losses of ice due to melt and sublimation)
   if(iSnow==1)then  ! if top snow layer include sublimation and melt
    volFracIceLoss = max(0._dp,mLayerMeltFreeze(iSnow)/iden_ice - dt*(scalarSnowSublimation/mLayerDepth(iSnow))/iden_ice )
   else
    volFracIceLoss = max(0._dp,mLayerMeltFreeze(iSnow)/iden_ice)  ! volumetric fraction of ice lost due to melt (-)
   end if
   ! (adjust snow depth to account for cavitation)
   scalarDepthNew = mLayerDepth(iSnow) * mLayerVolFracIceNew(iSnow)/(mLayerVolFracIceNew(iSnow) + volFracIceLoss)
   !print*, 'volFracIceLoss = ', volFracIceLoss
  else
   scalarDepthNew = mLayerDepth(iSnow)
  end if
  ! compute the total compaction rate associated with metamorphism
  CR_metamorph = CR_grainGrowth + CR_ovrvdnPress
  ! update depth due to metamorphism (implicit solution)
  mLayerDepth(iSnow) = scalarDepthNew/(1._dp + CR_metamorph*dt)

  ! check that depth is reasonable
  if(mLayerDepth(iSnow) < 0._dp)then
   write(*,'(a,1x,i4,1x,10(f12.5,1x))') 'iSnow, dt, density, massIceOld, massLiqOld = ', iSnow, dt, mLayerVolFracIceNew(iSnow)*iden_ice, massIceOld, massLiqOld
   write(*,'(a,1x,i4,1x,10(f12.5,1x))') 'iSnow, mLayerDepth(iSnow), scalarDepthNew, mLayerVolFracIceNew(iSnow), mLayerMeltFreeze(iSnow), CR_grainGrowth*dt, CR_ovrvdnPress*dt = ', &
                                         iSnow, mLayerDepth(iSnow), scalarDepthNew, mLayerVolFracIceNew(iSnow), mLayerMeltFreeze(iSnow), CR_grainGrowth*dt, CR_ovrvdnPress*dt
  endif

  ! update volumetric ice and liquid water content
  mLayerVolFracIceNew(iSnow) = massIceOld/(mLayerDepth(iSnow)*iden_ice)
  mLayerVolFracLiqNew(iSnow) = massLiqOld/(mLayerDepth(iSnow)*iden_water)

  !write(*,'(a,1x,i4,1x,f9.3)') 'after compact: iSnow, density = ', iSnow, mLayerVolFracIceNew(iSnow)*iden_ice
  !if(mLayerMeltFreeze(iSnow) > 20._dp) pause 'meaningful melt'

 end do  ! looping through snow layers

 ! check depth
 if(any(mLayerDepth(1:nSnow) < 0._dp))then
  do iSnow=1,nSnow
   write(*,'(a,1x,i4,1x,4(f12.5,1x))') 'iSnow, mLayerDepth(iSnow)', iSnow, mLayerDepth(iSnow)
  end do
  message=trim(message)//'unreasonable value for snow depth'
  err=20; return
 end if

 ! check for low/high snow density
 if(any(mLayerVolFracIceNew(1:nSnow)*iden_ice < minLayerDensity) .or. &
    any(mLayerVolFracIceNew(1:nSnow) > 1._dp))then
  do iSnow=1,nSnow
   write(*,'(a,1x,i4,1x,f9.3)') 'iSnow, density = ', iSnow, mLayerVolFracIceNew(iSnow)*iden_ice
  end do
  message=trim(message)//'unreasonable value for snow density'
  err=20; return
 end if

 end subroutine snwDensify


end module snwDensify_module
