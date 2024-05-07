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
! Model to calculate the dielectric constant ice
!   1. Hallikainen et al. 1995
!---------------------------------------------------------------------------

!SUBROUTINE DIEL_ICE 

! Purpose :
!   Calculate the relative permittivty of pure ice in the microwave region

! Reference:
!  P60-68, ESTEC CONTRACT, NO 11706/95/NL/NB(SC), FINAL REPORT, HALLIKAINEN ET AL.
! CODE DOWNLOADED FROM HUT SITE (J. POULLIAINEN)  - 2001
!       EPX=3.15,  EPY=-1e-3,  EPSIC=CMPLX(EPX,EPY) (ULABY ET AL., 1986, P2026)    

! Input/Output Variables:
! T (K)
! eps_ice : dielectric constant of ice
! Input Variables:
!   fghz         frequency (GHz)
!   T            temperature (K) 

! Output variables:
!   eps_ice      dielectric constant of ice 


! Internal variables:
!   eir		 real part of dielectric constant of ice (MATZLER & WEGMULLER 1987)
!   tfreeze      (K) 
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------

SUBROUTINE DIEL_ICE(fghz, T, eps_ice) 

IMPLICIT NONE

real :: fghz, T
complex :: eps_ice

REAL :: B1, B2, BB, DELTABETA, BETAM, BETA,THETA,ALFA
REAL :: eir, tfreeze
COMPLEX :: E
!---------------------------------------------------------------------------

tfreeze = 273.15
E=(0.,1.)

! -------- IMPROVED VERSION (MISHIMA, MATZLER)
B1 = 0.0207
B2 = 1.16e-11
BB = 335.
DELTABETA = EXP(-10.02 + 0.0364*(T-tfreeze))
BETAM = (B1/T) * ( EXP(BB/T) / ((EXP(BB/T)-1.)**2.) ) + B2*(fghz**2.)
BETA = BETAM + DELTABETA
! -------- PURE ICE AS A FUNCTION OF TEMPERATURE (HUFFORD 1991):
THETA = 300. / T - 1.
ALFA = (0.00504 + 0.0062 * THETA) * EXP(-22.1 * THETA)

! -------- REAL PART OF ICE EPSILON (MATZLER & WEGMULLER 1987)
eir = 3.1884 + 9.1e-4 * (T-tfreeze)

eps_ice = eir - E*(ALFA/fghz + BETA *fghz)
   
END SUBROUTINE DIEL_ICE
