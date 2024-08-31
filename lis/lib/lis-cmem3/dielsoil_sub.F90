!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!!!! Insert banner here !!!!
!---------------------------------------------------------------------------
! Models to calculate the dielectric constant of a soil medium
!   1. Wang and Schmugge, 1980
!   2. Dobson et al., 1985:
!---------------------------------------------------------------------------

!SUBROUTINE DIELWANG 

! Purpose :
!   Calculate the dielectric constant of a wet soil 
!   Developed and validated for 1.4 and 5 GHz.

! Reference:
!  Wang and Schmugge, 1980: An empirical model for the 
!    complex dielectric permittivity of soils as a function of water
!    content. IEEE Trans. Geosci. Rem. Sens., GE-18, No. 4, 288-295.

! External :
!   dielsal_sub.F90

! Input Variables:
!   fghz	 frequency (GHz) 
!   ew		 dielectric constant of water   
!   wc		 volumetric soil water content (cm3/cm3)
!   sand  	 sand fraction (% of dry weight of soil)
!   clay  	 clay fraction (% of dry weight of soil)
!   p 		 porosity of dry soil (cm3/cm3)

! Output variables: 
!   eps		 dielectric constant of the wet soil  


! Internal variables: 
!   ei 		 dielectric constant of ice, for initially absorbed water
!   ea  	 dielectric constant of air
!   er  	 dielectric constant of rock 
!   wp		 wilting point (cm3/cm3) 
!   alpha 	 fitting parameter (Wang and Schmugge, 1980)
!   wt : transition moisture point (cm3/cm3)
!   gamma : fitting parameter
!   ecl : conductivity loss
!   ex : dielectric constant of the initially absorbed water
!---------------------------------------------------------------------------

SUBROUTINE DIELWANG(fghz, ew, wc, sand, clay, p, eps) 

IMPLICIT NONE

real :: fghz, wc, sand, clay, p
COMPLEX :: ew, eps
REAL :: wt, wp, alpha, gamma, ecl
COMPLEX ::  ei, ea, er, ex 
!---------------------------------------------------------------------------

ei = (3.2, 0.1)
ea = (1.0, 0.0) 
er = (5.5, 0.2)
wp = 0.06774 - 0.00064 * sand + 0.00478 * clay

if ( fghz > 2.5 ) then
        alpha=0.
else
        alpha=100. * wp 
        if ( alpha > 26.) alpha=26.
endif 

! 1. Calculate dielectric constant of soil-water mixture

gamma = -0.57 * wp + 0.481
wt = 0.49 * wp + 0.165
IF (wc <= wt) THEN
  ex = ei + (ew-ei) * (wc/wt) * gamma
  eps = wc*ex + (p-wc) * ea + (1.-p) * er
ELSE
  ex = ei + (ew-ei) * gamma
  eps = wt*ex + (wc-wt) * ew + (p-wc) * ea + (1.-p) * er
ENDIF

! 2. add conductivity loss

ecl = alpha * wc**2.
eps = eps + (0.,1.) * ecl

END SUBROUTINE DIELWANG

!===========================================================================

!SUBROUTINE DIELDOBSON 

! Purpose :
!   Calculate the dielectric constant of a wet soil 
!   Developed and validated for 1.4 and 18 GHz.

! Reference:
!  Dobson et al., 1985: Microwave Dielectric behavior of 
!    wet soil - part II: Dielectric mixing models,
!    IEEE Trans. Geosc. Rem. Sens., GE-23, No. 1, 35-46.


! Input Variables:
!   ew           dielectric constant of water
!   wc           volumetric soil water content (cm3/cm3)
!   sand         sand fraction (% of dry weight of soil)
!   clay         clay fraction (% of dry weight of soil)

! Output variables:
!   eps          dielectric constant of the wet soil


! Internal variables:
!   rho_s	 rho_s : soil specific density (g/cm3), SMOS ATBD 2.664
!   rho_b	 soil bulk density (g/cm3) ~1.3
!   alphas 	 constant shape factor
!---------------------------------------------------------------------------

SUBROUTINE DIELDOBSON (ew, wc, sand, clay, eps)

IMPLICIT NONE

REAL :: wc, sand, clay
COMPLEX :: ew, eps
REAL :: rho_s, rho_b, eps_s, beta, eaa,  epsr, epsi
REAL, PARAMETER :: alphas = 0.65
!---------------------------------------------------------------------------

rho_s = 2.66
rho_b = (sand*1.6 +  clay*1.1 + (100.-sand-clay)*1.2)/100.

wc = MAX(wc,0.001) ! to avoid dividing by zero

eps_s = (1.01 + 0.44 * rho_s)**2. - 0.062  ! compare to ATBD: 4.7

beta = (127.48 - 0.519 * sand - 0.152 * clay) / 100.
eaa = 1.0 + (rho_b / rho_s) * (eps_s ** alphas - 1.0)   &
  & + (wc ** beta) * (real(ew) ** alphas) - wc
epsr = eaa ** (1./alphas)

beta = (133.797 - 0.603 * sand - 0.166 * clay) / 100.
eaa= (wc ** beta) * (abs(aimag(ew)) ** alphas)
epsi = eaa ** (1./alphas)

eps = cmplx(epsr,epsi)

END SUBROUTINE DIELDOBSON


!=====================================================================================

!SUBROUTINE DIELMIRONOV 

! Purpose :
!   Calculate the dielectric constant of a wet soil
!   Developed and validated from 1 to 10 GHz.
!   adapted for a large range of soil moisture

! Reference:
!   Mironov et al: Generalized Refractive Mixing Dielectric Model for moist soil
!    IEEE Trans. Geosc. Rem. Sens., vol 42 (4), 773-785. 2004.
!
! Adapted from the Matlab version of JP Wigneron
! 
! Patricia de Rosnay, October 9 2007

! Input Variables:
!   fghz	 frequency (GHz) 
!   wc           volumetric soil water content (cm3/cm3)
!   clay         clay fraction (% of dry weight of soil)

! Output variables:
!   eps          dielectric constant of the wet soil


! Internal variables:
!   f		 frequency (Hz) 
!   eps_winf 	 dielectric constant at infinite frequency (Stogryn 1971)
!       	 also called: high frequency dielectric constant (4.9 Lane and Saxton )
!   eps_0 	 permittivity of free space (Klein and Swift 1977) [Farads/meter]
!---------------------------------------------------------------------------

SUBROUTINE DIELMIRONOV (fghz, wc, clay, eps) 

IMPLICIT NONE

real :: fghz, wc, clay
complex :: eps

real :: f, pi, eps_0, eps_winf
REAL :: znd,zkd,zxmvt,zep0b,ztaub,zsigmab,zep0u,ztauu,zsigmau,zcxb,zepwux,zepwuy &
                    &, zcxu,zepwbx,zepwby,znb,zkb,znu,zku,znm,zkm,zflag,zxmvt2,zepmx,zepmy 

f = fghz * 1.0e9
pi = acos(-1.0) 
eps_0 = 8.854e-12
eps_winf = 4.9

! Initializing the GRMDM spectroscopic parameters with clay (fraction)
!-----------------------------------------------------------------------

!  RI & NAC of dry soils
!-----------------------
znd = 1.634 - 0.539 * (clay/100.) + 0.2748 * (clay/100.)**2.
zkd = 0.03952 - 0.04038 * (clay / 100.)

! Maximum bound water fraction
!-----------------------------
zxmvt = 0.02863 + 0.30673 * clay / 100.

! Bound water parameters
!-----------------------
zep0b = 79.8 - 85.4  * (clay / 100.) + 32.7  * (clay / 100.)*(clay / 100.)
ztaub = 1.062e-11 + 3.450e-12 * (clay / 100.)
zsigmab = 0.3112 + 0.467 * (clay / 100.)

! Unbound (free) water parameters
!-------------------------------
zep0u = 100.
ztauu = 8.5e-12
zsigmau = 0.3631 + 1.217 * (clay / 100.)


! Computation of epsilon water (bound & unbound)
!----------------------------------------------

zcxb = (zep0b - eps_winf) / (1. + (2.*pi*f*ztaub)**2.)
zepwbx = eps_winf + zcxb
zepwby =  zcxb * (2.*pi*f*ztaub) + zsigmab / (2.*pi*eps_0*f)

zcxu = (zep0u - eps_winf) / (1. + (2.*pi*f*ztauu)**2.)
zepwux = eps_winf + zcxu
zepwuy =  zcxu * (2.*pi*f*ztauu) + zsigmau/(2.*pi*eps_0*f)

! Computation of refractive index of water (bound & unbound)
! -----------------------------------------------------------

znb= sqrt( sqrt( zepwbx**2. + zepwby**2.) + zepwbx ) / sqrt(2.)
zkb= sqrt( sqrt( zepwbx**2. + zepwby**2.) - zepwbx ) / sqrt(2.)

znu= sqrt( sqrt( zepwux**2. + zepwuy**2.) + zepwux ) / sqrt(2.)
zku= sqrt( sqrt( zepwux**2. + zepwuy**2.) - zepwux ) / sqrt(2.)

! Computation of soil refractive index (nm & km):
! xmv can be a vector
! -----------------------------------------------------------

zxmvt2= min (wc, zxmvt)
zflag = 0.
IF ( wc >= zxmvt ) zflag = 1.

znm = znd + (znb - 1.) * zxmvt2 + (znu - 1.) * (wc-zxmvt) * zflag
zkm = zkd + zkb * zxmvt2 + zku * (wc-zxmvt) * zflag

! computation of soil dielectric constant:
! -----------------------------------------------------------

zepmx= znm ** 2. - zkm ** 2.     
zepmy= znm * zkm * 2.      

eps = cmplx(zepmx,zepmy)


END SUBROUTINE DIELMIRONOV


