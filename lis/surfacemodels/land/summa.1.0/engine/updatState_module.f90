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

module updatState_module
USE nrtype
! physical constants
USE multiconst,only:&
                    Tfreeze,     & ! freezing point of pure water  (K)
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  & ! intrinsic density of water    (kg m-3)
                    gravity,     & ! gravitational acceleteration  (m s-2)
                    LH_fus         ! latent heat of fusion         (J kg-1)
implicit none
private
public::updateSnow
public::updateSoil
contains


 ! *************************************************************************************************************
 ! public subroutine updateSnow: compute phase change impacts on volumetric liquid water and ice
 ! *************************************************************************************************************
 subroutine updateSnow(&
                       ! input
                       mLayerTemp       ,& ! intent(in): temperature (K)
                       mLayerTheta      ,& ! intent(in): volume fraction of total water (-)
                       snowfrz_scale    ,& ! intent(in): scaling parameter for the snow freezing curve (K-1)
                       ! output
                       mLayerVolFracLiq ,& ! intent(out): volumetric fraction of liquid water (-)
                       mLayerVolFracIce ,& ! intent(out): volumetric fraction of ice (-)
                       fLiq             ,& ! intent(out): fraction of liquid water (-)
                       err,message)        ! intent(out): error control
 ! utility routines
 USE snow_utils_module,only:fracliquid     ! compute volumetric fraction of liquid water
 implicit none
 ! input variables
 real(dp),intent(in)           :: mLayerTemp           ! temperature (K)
 real(dp),intent(in)           :: mLayerTheta          ! volume fraction of total water (-)
 real(dp),intent(in)           :: snowfrz_scale        ! scaling parameter for the snow freezing curve (K-1)
 ! output variables
 real(dp),intent(out)          :: mLayerVolFracLiq     ! volumetric fraction of liquid water (-)
 real(dp),intent(out)          :: mLayerVolFracIce     ! volumetric fraction of ice (-)
 real(dp),intent(out)          :: fLiq                 ! fraction of liquid water (-)
 ! error control
 integer(i4b),intent(out)      :: err                  ! error code
 character(*),intent(out)      :: message              ! error message
 ! initialize error control
 err=0; message="updateSnow/"

 ! compute the volumetric fraction of liquid water and ice (-)
 fLiq = fracliquid(mLayerTemp,snowfrz_scale)
 mLayerVolFracLiq = fLiq*mLayerTheta
 mLayerVolFracIce = (1._dp - fLiq)*mLayerTheta*(iden_water/iden_ice)
 !print*, 'mLayerTheta - (mLayerVolFracIce*(iden_ice/iden_water) + mLayerVolFracLiq) = ', mLayerTheta - (mLayerVolFracIce*(iden_ice/iden_water) + mLayerVolFracLiq)

 !write(*,'(a,1x,4(f20.10,1x))') 'in updateSnow: fLiq, mLayerTheta, mLayerVolFracIce = ', &
 !                                               fLiq, mLayerTheta, mLayerVolFracIce
 !pause

 end subroutine updateSnow

 ! *************************************************************************************************************
 ! public subroutine updateSoil: compute phase change impacts on matric head and volumetric liquid water and ice
 ! *************************************************************************************************************
 subroutine updateSoil(&
                       ! input
                       mLayerTemp       ,& ! intent(in): temperature vector (K)
                       mLayerMatricHead ,& ! intent(in): matric head (m)
                       vGn_alpha        ,& ! intent(in): van Genutchen "alpha" parameter
                       vGn_n            ,& ! intent(in): van Genutchen "n" parameter
                       theta_sat        ,& ! intent(in): soil porosity (-)
                       theta_res        ,& ! intent(in): soil residual volumetric water content (-)
                       vGn_m            ,& ! intent(in): van Genutchen "m" parameter (-)
                       ! output
                       mLayerVolFracWat ,& ! intent(out): volumetric fraction of total water (-)
                       mLayerVolFracLiq ,& ! intent(out): volumetric fraction of liquid water (-)
                       mLayerVolFracIce ,& ! intent(out): volumetric fraction of ice (-)
                       err,message)        ! intent(out): error control
 ! utility routines
 USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water based on matric head
 USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric liquid water content
 implicit none
 ! input variables
 real(dp),intent(in)           :: mLayerTemp           ! estimate of temperature (K)
 real(dp),intent(in)           :: mLayerMatricHead     ! matric head (m)
 real(dp),intent(in)           :: vGn_alpha            ! van Genutchen "alpha" parameter
 real(dp),intent(in)           :: vGn_n                ! van Genutchen "n" parameter
 real(dp),intent(in)           :: theta_sat            ! soil porosity (-)
 real(dp),intent(in)           :: theta_res            ! soil residual volumetric water content (-)
 real(dp),intent(in)           :: vGn_m                ! van Genutchen "m" parameter (-)
 ! output variables
 real(dp),intent(out)          :: mLayerVolFracWat     ! fractional volume of total water (-)
 real(dp),intent(out)          :: mLayerVolFracLiq     ! volumetric fraction of liquid water (-)
 real(dp),intent(out)          :: mLayerVolFracIce     ! volumetric fraction of ice (-)
 integer(i4b),intent(out)      :: err                  ! error code
 character(*),intent(out)      :: message              ! error message
 ! define local variables
 real(dp)                      :: TcSoil               ! critical soil temperature when all water is unfrozen (K)
 real(dp)                      :: xConst               ! constant in the freezing curve function (m K-1)
 real(dp)                      :: mLayerPsiLiq         ! liquid water matric potential (m)
 ! initialize error control
 err=0; message="updateSoil/"

 ! compute fractional **volume** of total water (liquid plus ice)
 mLayerVolFracWat = volFracLiq(mLayerMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
 if(mLayerVolFracWat > theta_sat)then; err=20; message=trim(message)//'volume of liquid and ice exceeds porosity'; return; end if

 ! compute the critical soil temperature where all water is unfrozen (K)
 ! (eq 17 in Dall'Amico 2011)
 TcSoil = Tfreeze + min(mLayerMatricHead,0._dp)*gravity*Tfreeze/LH_fus  ! (NOTE: J = kg m2 s-2, so LH_fus is in units of m2 s-2)

 ! *** compute volumetric fraction of liquid water and ice for partially frozen soil
 if(mLayerTemp < TcSoil)then ! (check if soil temperature is less than the critical temperature)

  ! - volumetric liquid water content (-)
  ! NOTE: mLayerPsiLiq is the liquid water matric potential from the Clapeyron equation, used to separate the total water into liquid water and ice
  !       mLayerPsiLiq is DIFFERENT from the liquid water matric potential used in the flux calculations
  xConst           = LH_fus/(gravity*Tfreeze)        ! m K-1 (NOTE: J = kg m2 s-2)
  mLayerPsiLiq     = xConst*(mLayerTemp - Tfreeze)   ! liquid water matric potential from the Clapeyron eqution
  mLayerVolFracLiq = volFracLiq(mLayerPsiLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)

  ! - volumetric ice content (-)
  mLayerVolFracIce = mLayerVolFracWat - mLayerVolFracLiq

 ! *** compute volumetric fraction of liquid water and ice for unfrozen soil
 else

  ! all water is unfrozen
  mLayerPsiLiq     = mLayerMatricHead
  mLayerVolFracLiq = mLayerVolFracWat
  mLayerVolFracIce = 0._dp

 end if  ! (check if soil is partially frozen)

 end subroutine updateSoil


end module updatState_module
