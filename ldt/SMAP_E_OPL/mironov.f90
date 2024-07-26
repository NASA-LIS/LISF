!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
MODULE mironov_m

        IMPLICIT NONE

      CONTAINS

        SUBROUTINE mironov(freq,mv,clayfrac,er_r)

          IMPLICIT NONE

          REAL(4), INTENT(IN)     :: freq, mv, clayfrac
          COMPLEX(4), INTENT(OUT) :: er_r

          REAL(4), PARAMETER      :: PI = acos(-1.0)
          REAL(4)                 :: f, C
          REAL(4)                 :: nd, kd, mvt, eps0b, eps0u, taub, tauu, sigb, sigu
          REAL(4)                 :: eps0, epsinf, epsb_real, epsu_real, epsb_imag, epsu_imag
          REAL(4)                 :: nb, kb, nu, ku
          REAL(4)                 :: nm, km
          REAL(4)                 :: er_r_real, er_r_imag

          f = freq*1e9                                                                              ! Section IV
          C = clayfrac*100                                                                          ! Section VI

       !! Mironov's regression expressions based on Curtis, Dobson, and Hallikainen datasets
          nd = 1.634 - 0.539e-2 * C + 0.2748e-4 * C**2                                              ! Eqn 17
          kd = 0.03952 - 0.04038e-2 * C                                                             ! Eqn 18
          mvt = 0.02863 + 0.30673e-2 * C                                                            ! Eqn 19
          eps0b = 79.8 - 85.4e-2 * C + 32.7e-4 * C**2                                               ! Eqn 20
          taub = 1.062e-11 + 3.450e-12 * 1e-2 * C                                                   ! Eqn 21
          sigb = 0.3112 + 0.467e-2 * C                                                              ! Eqn 22
          sigu = 0.3631 + 1.217e-2 * C                                                              ! Eqn 23
          eps0u = 100                                                                               ! Eqn 24
          tauu = 8.5e-12                                                                            ! Eqn 25

       !! Debye relaxation equations for water as a function of frequency
          eps0 = 8.854e-12                                                                          ! Vacuum permittivity
          epsinf = 4.9                                                                              ! Section IV
          epsb_real = epsinf + ( (eps0b - epsinf)/(1 + (2*PI*f*taub)**2) )                          ! Eqn 16
          epsb_imag = (eps0b - epsinf)/(1 + (2*PI*f*taub)**2) * (2*PI*f*taub) + sigb/(2*PI*eps0*f)  ! Eqn 16
          epsu_real = epsinf + ( (eps0u - epsinf)/(1 + (2*PI*f*tauu)**2) )                          ! Eqn 16
          epsu_imag = (eps0u - epsinf)/(1 + (2*PI*f*tauu)**2) * (2*PI*f*tauu) + sigu/(2*PI*eps0*f)  ! Eqn 16

       !! Refractive indices
          nb = 1/sqrt(2.0) * sqrt( sqrt(epsb_real**2 + epsb_imag**2) + epsb_real )                  ! Eqn 14
          kb = 1/sqrt(2.0) * sqrt( sqrt(epsb_real**2 + epsb_imag**2) - epsb_real )                  ! Eqn 14
          nu = 1/sqrt(2.0) * sqrt( sqrt(epsu_real**2 + epsu_imag**2) + epsu_real )                  ! Eqn 14
          ku = 1/sqrt(2.0) * sqrt( sqrt(epsu_real**2 + epsu_imag**2) - epsu_real )                  ! Eqn 14

       !! n(*) are refractive indices, k(*) are normalized attenuation coefficients
       !!   m: moist soil
       !!   d: dry soil
       !!   b: bound soil water (BSW)
       !!   u: unbound (free) soil water (FSW)
          IF (mv <= mvt) THEN
              nm = nd + (nb - 1) * mv                                                               ! Eqn 12
              km = kd + kb * mv                                                                     ! Eqn 13
              er_r_real = nm**2 - km**2                                                             ! Eqn 11
              er_r_imag = 2 * nm * km                                                               ! Eqn 11
          ELSE
              nm = nd + (nb - 1) * mvt + (nu - 1) * (mv - mvt)                                      ! Eqn 12
              km = kd + kb * mvt + ku * (mv - mvt)                                                  ! Eqn 13
              er_r_real = nm**2 - km**2                                                             ! Eqn 11
              er_r_imag = 2 * nm * km                                                               ! Eqn 11
          ENDIF

          er_r = CMPLX(er_r_real,er_r_imag)
        END SUBROUTINE mironov

      END MODULE mironov_m

