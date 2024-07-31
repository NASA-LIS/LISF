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
! Different models to calculate Rough Surface Emissivity
!  1,2,4,5) Choudhury et al., 1979
!  3) Wegmueller and Maetzler, 1999
!---------------------------------------------------------------------------

!SUBROUTINE RGHCHOU

! Purpose :
!   Calculate Rough Surface Emissivity

! Reference:
!    Wang and Choudhury, 1981: Remote sensing of soil moisture
!     content over bare field at 1.4 GHz frequency, J.Geo.Res.
!     Vol.86, 5277-5287    
!    Choudhury et al., 1979: Effect of surface roughness on the
!     microwave emission from soils, J.Geo.Res. Vol.84, 5699-5706
! and what reference for the parameterization of h and Q?

! Input variables: 
! fghz          frequency (GHz) 
! theta         incidence angle (degrees)
! Nrh		exponent in surface roughness expression, h-pol 
!               function of land cover 
! Nrv 		exponent in surface roughness expression, v-pol
!               function of land cover 
! hrmodel     	vegitation-determined parameter 
!               function of land cover 
! r_s		smooth surface reflectivity 

! Output variables: 
! r_r		rought surface reflectivity 
!---------------------------------------------------------------------------

SUBROUTINE RGHCHOU(fghz, theta, Nrh, Nrv, hrmodel, r_s, r_r)

IMPLICIT NONE
real :: fghz, theta, Nrh, Nrv, hrmodel, r_s(2), r_r(2)

real :: pi, costheta, ip_rgh_surf, ip_Q

!---------------------------------------------------------------------------
pi = acos(-1.0)
costheta = cos(theta*pi/180.0)

! ip_Q : Cross polarization parameter:

! CMEM original, CRTM : 0.5
!ip_rgh_surf = 2.2
! CRTM
ip_rgh_surf = 0.5 


ip_Q = 0.35 * (1.0 - exp(-0.6 * ip_rgh_surf**2 * fghz))
IF ( fghz < 2. ) ip_Q=0.   ! Q is assumed zero at low frequency


! CMEM original
!r_r(1) = (ip_Q * r_s(2) + (1.-ip_Q) * r_s(1)) * exp(-hrmodel * costheta**Nrh)
!r_r(2) = (ip_Q * r_s(1) + (1.-ip_Q) * r_s(2)) * exp(-hrmodel * costheta**Nrv)

! CRTM formulation
r_r(1) = 0.3 * (ip_Q * r_s(2) + (1.-ip_Q) * r_s(1)) 
r_r(2) = 0.3 * (ip_Q * r_s(1) + (1.-ip_Q) * r_s(2)) 

END SUBROUTINE RGHCHOU


! Same as RGHCHOU, but use special parameters for desert

SUBROUTINE RGHdesert(fghz, theta, Nrh, Nrv, hrmodel, r_s, r_r)

IMPLICIT NONE
real :: fghz, theta, Nrh, Nrv, hrmodel, r_s(2), r_r(2)

real :: pi, costheta, ip_rgh_surf, ip_Q

!---------------------------------------------------------------------------
pi = acos(-1.0)
costheta = cos(theta*pi/180.0)

! ip_Q : Cross polarization parameter:

! CMEM original, CRTM : 0.5
!ip_rgh_surf = 2.2
! CRTM
ip_rgh_surf = 0.1


ip_Q = 0.35 * (1.0 - exp(-0.6 * ip_rgh_surf**2 * fghz))
IF ( fghz < 2. ) ip_Q=0.   ! Q is assumed zero at low frequency


! CMEM original
!r_r(1) = (ip_Q * r_s(2) + (1.-ip_Q) * r_s(1)) * exp(-hrmodel * costheta**Nrh)
!r_r(2) = (ip_Q * r_s(1) + (1.-ip_Q) * r_s(2)) * exp(-hrmodel * costheta**Nrv)

! CRTM formulation
r_r(1) = 0.3 * (ip_Q * r_s(2) + (1.-ip_Q) * r_s(1))
r_r(2) = 0.3 * (ip_Q * r_s(1) + (1.-ip_Q) * r_s(2))

END SUBROUTINE RGHdesert


!===========================================================================

!SUBROUTINE RGHWEGM

! Purpose :
! calculate Rough Surface Emissivity based on smooth surface reflectivity (H only)

!  Reference:
!    Wegmueller and Maetzler, 1999: Rough bare soil reflectivity model,
!    IEEE Trans. Geosci. Remote Sensing, Vol.37, No.3, 1391-1395.

! Input variables:
! fghz          frequency (Ghz) 
! theta         incidence angle (degrees)
! s  		effective roughness height [m], empirical function of land cover
! r_s           smooth surface reflectivity

! Output variables:
! r_r           rought surface reflectivity


!---------------------------------------------------------------------------
SUBROUTINE RGHWEGM(fghz, theta, s, r_s, r_r)

IMPLICIT NONE
real :: fghz, theta, s, r_s(2), r_r(2) 
real :: pi, costheta, k, ks
real, parameter :: c = 2.998E8        ! speed of light m/s

!---------------------------------------------------------------------------
pi = acos(-1.0)
costheta = cos(theta*pi/180.0)

k = 2.0 * pi * fghz * 1.0E9 / c         !  wave number, 1/m
ks = k * s                              !  1/m * m 

r_r(1) = r_s(1) * exp(-1.0 * ks ** (sqrt(0.10 * costheta)))

IF (ks == 0.) THEN
  r_r(2) = r_s(2) 
ELSE
  IF (theta <= 60.0) THEN
    r_r(2) = r_r(1) * (costheta ** 0.655)
  ELSEIF (theta <= 70.0) THEN
    r_r(2) = r_r(1) * (0.635 - 0.0014*(theta-60.0))
  ELSE
    Write(*, *)'Incidence angle not in [0;70] degrees. not suitable for rough soil reflectivity'
    Stop
  ENDIF
ENDIF

END SUBROUTINE RGHWEGM
