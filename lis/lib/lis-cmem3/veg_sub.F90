!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! vegetation opacity models (iveg)
!               0 = no vegetation
!               1 = Kirdyashev et al., 1979
!               2 = Wegmueller et al., 1995
!               3 = Wigneron
!---------------------------------------------------------------------------
!  About vegetation structure coefficient a_geo:
!
!    a_geo(pol) = u(pol) / 3
!
!    u takes vegetation geometry into account
!
!    Theoretical examples:
!       horizontal leaves:  a_geo(h) = 1     a_geo(v) = (cos(theta))**2
!       isotropic leaves:   a_geo(h) = 2/3   a_geo(v) = 2/3
!       vertical stalks:    a_geo(h) = 0     a_geo(v) = (sin(theta))**2
!       isotropic stalks:   a_geo(h) = 1/3   a_geo(v) = 1/3
!---------------------------------------------------------------------------

!SUBROUTINE VEGWEGM

! Purpose :
!  Calculate vegetation opacity using 'Geometrical Optics theory'
!  For low to high  frequencies (<100GHz)

!  Reference:
!  Wegmueller et al., 1995: Canopy opacity models, Passive Microwave remote 
!    sensing of land-atmosphere interactions, VSP, 375 ff.
! 1. Matzler, ??

! Input Variables:
!   fghz         frequency (GHz)
!   theta        zenith angle (degree) 
!   eps_vw       dielectric constant of vegetation water
!   wc_veg       vegetation water content per unit area (kg/m2)
!   a_geo	 vegetation structure coefficient (index 1.h-pol, 2.v-pol)
!   rho_veg 	 vegetation density (kg/m3)
!   m_d     	 dry mass fraction of vegetation
!   d_leaf  	 leaf thickness (mm)

! Output Variables:
!   tau_veg      effective vegetation opacity (index 1.h-pol, 2.v-pol)
!   omega_veg    single scattering albedo (index 1.h-pol, 2.v-pol)

! Internal Variables:
!  eps_veg : Dielectric constant of leaves
!  BB      : wet biomass
!---------------------------------------------------------------------------

SUBROUTINE VEGWEGM(fghz, theta, &                ! sensor params
                   eps_vw, wc_veg, a_geo,   &    ! veg params
                   rho_veg, m_d, d_leafmm,  & 
                   tau_veg, omega_veg)           ! output
IMPLICIT NONE

real :: fghz, theta
real :: wc_veg, a_geo(2), rho_veg, m_d, d_leafmm, d_leaf
real :: tau_veg(2), omega_veg(2)
complex :: eps_vw

! wave parameters 
real, parameter :: c = 2.998E8            ! light speed,  m/s
real :: pi, f, k, costheta, sintheta

COMPLEX :: eps_veg
REAL :: kz0
REAL :: th, tv, rh, rv
REAL :: BB
COMPLEX :: kz1
COMPLEX :: RRh, RRv, j, hfrac, vfrac, exp1, exp2

 pi = acos(-1.0)
 f = fghz * 1.0E9
 k = 2.0 * pi * f / c
 costheta = cos(theta*pi/180.0)
 sintheta = sin(theta*pi/180.0)

! CMEM default values 
! rho_veg = 950.
! m_d = 0.3
! d_leafmm = 0.2

d_leaf = d_leafmm * 0.001     ! convert mm to m
!---------------------------------------------------------------------------

eps_veg = (0.51 - 0.64 * m_d) * eps_vw  +  3.2 * m_d  +  0.49

! 1. Calculate transmissivity of a single leaf
!    -----------------------------------------
kz0 = k * costheta
kz1 = k * sqrt(eps_veg - sintheta * sintheta)

RRh = (kz0 - kz1) / (kz0 + kz1)
RRv = (eps_veg * kz0 - kz1) / (eps_veg * kz0 + kz1)

j = (0., -1.)    ! YDT: this is CMEM original. Can not use j = (0., 1.) -- it will produce
                 !      insane values. 
exp1 = exp(j * (kz0 - kz1) * d_leaf)
exp2 = exp(-j * 2. * kz1 * d_leaf)

hfrac = RRh * ( 1.0 - exp2)/ ( 1.0 - RRh*RRh * exp2 )
vfrac = RRv * ( 1.0 - exp2)/ ( 1.0 - RRv*RRv * exp2 )

#if (defined DEBUG)
write(*, *) "veg_sub.F90: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"  abs(RRh)  abs(RRv) abs(hfrac) abs(vfrac) abs(exp1) abs(exp2)" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(6F10.3)') abs(RRh), abs(RRv), abs(hfrac), abs(vfrac), abs(exp1), abs(exp2) 
write(*, *)
#endif

rh = abs(hfrac) * abs(hfrac) 
rv = abs(vfrac) * abs(vfrac) 

hfrac = (4. * kz0 * kz1 * exp1) / ((kz0 + kz1)**2 * (1. - RRh*RRh * exp2))
vfrac = (4. * eps_veg * kz0 * kz1 * exp1) /   &
        ((eps_veg * kz0 + kz1)**2 * (1. - RRv*RRv * exp2))

th = abs(hfrac) * abs(hfrac) 
tv = abs(vfrac) * abs(vfrac) 

! 1. Calculate single scattering albedo 
!    ----------------------------
omega_veg(1) = rh / ( 1.0 - th) 
omega_veg(2) = rv / ( 1.0 - tv) 

#if (defined DEBUG)
write(*, *) "veg_sub.F90: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"                fghz     theta        rh        rv        th        tv" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(A10, 6F10.3)') "r:t", fghz, theta, rh, rv, th, tv
write(*, *)
#endif

a_geo(1) = 1.0
a_geo(2) = 1.0 
! 2. Calculate vegetation opacity
!    ----------------------------
BB = wc_veg / (1. - m_d)   

tau_veg(1) = a_geo(1) * k * ( BB / rho_veg ) * aimag(eps_veg) * th / costheta
tau_veg(2) = a_geo(2) * k * ( BB / rho_veg )  * aimag(eps_veg) * tv / costheta

! CRTM approach
tau_veg(1) = tau_veg(2) 

#if (defined DEBUG)
write(*, *) "veg_sub.F90: "
!123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
write(*, *)  &
"                fghz     theta    omegaH    omegaV  tau_vegH  tau_vegV"
write(*, *)"---------------------------------------------------------------------------"
write(*, '(A10, 6F10.3)') "omg:tau", fghz, theta, omega_veg(1), omega_veg(2), &
                          tau_veg(1), tau_veg(2)
write(*, *)  &
"   eps_veg              eps_vw" 
write(*, *)"---------------------------------------------------------------------------"
write(*, '(2(A, F7.1, A, F7.1, A))') & 
          '(', real(eps_veg) , ' + i', imag(eps_veg),  ') ',  &
          ' (', real(eps_vw) , ' + i', imag(eps_vw),  ') ' 
write(*, *)
#endif

!WRITE(6,*) wc_veg, m_d, eps_vw, eps_veg , ( B / rho_veg )

END SUBROUTINE VEGWEGM

