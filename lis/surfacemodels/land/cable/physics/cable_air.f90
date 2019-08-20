! cable_air.f90
!
! Source file containing air routines for CABLE
!
! Development by Ying-Ping Wang, Eva Kowalczyk, Ray Leuning
! Gab Abramowitz, Martin Dix, Harvey Davies, Mike Raupach, 
!
! bugs to bernard.pak@csiro.au
!
! This file contains modules:
!   air_module, 
! The subroutines included are:
!   define_air,
!
! Most user-defined types (e.g. met%tk) are defined in cable_types module
! in cable_variables.f90

MODULE cable_air
  USE cable_physical_constants
  USE cable_dimensions, ONLY: mp_patch, r_1
  USE cable_types
  IMPLICIT NONE
  PRIVATE
  PUBLIC define_air
CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE define_air(met,air)
    TYPE (air_type), INTENT(INOUT) :: air ! air_type variables
    TYPE (met_type), INTENT(IN)  :: met ! meteorological variables
    REAL(r_1), DIMENSION(mp_patch) :: es ! sat vapour pressure (mb)    
    ! Calculate saturation vapour pressure
    es = tetena * EXP(tetenb * (met%tvair-tfrz)/(tetenc + (met%tvair-tfrz)))
    ! Calculate conversion factor from from m/s to mol/m2/s
    air%cmolar = met%pmb * 100.0 / (rgas * (met%tvair))
    ! Calculate dry air density:
    air%rho = MIN(1.3,rmair * air%cmolar)
    ! molar volume (m^3/mol)
    air%volm = rgas * (met%tvair) / (100.0 * met%pmb)
    ! latent heat for water (j/kg)
    ! air%rlam = (2501.0 - 2.38 * (met%tvair- tfrz)) * 1000.0
	air%rlam  = 1.91846e6*((met%tvair)/((met%tvair)-33.91))**2  ! Henderson-Sellers, QJRMS, 1984
    !air%rlam=2.5104e6
    ! saturation specific humidity
    air%qsat = (rmh2o / rmair) * es / met%pmb
    ! d(qsat)/dT ((kg/kg)/K)
    air%epsi = (air%rlam / capp) * (rmh2o / rmair) * es * tetenb * tetenc / &
         &  (tetenc + (met%tvair - tfrz)) ** 2 / met%pmb
    ! air kinematic viscosity (m^2/s)
    air%visc = 1e-5 * MAX(1.0, 1.35 + 0.0092 * (met%tvair - tfrz))
    ! psychrometric constant
    air%psyc = met%pmb * 100.0 * capp * rmair / air%rlam / rmh2o
    ! d(es)/dT (mb/K)
    air%dsatdk = 100.0*(tetena*tetenb*tetenc)/((met%tvair-tfrz)+tetenc)**2 &
         * EXP(tetenb*(met%tvair-tfrz)/((met%tvair-tfrz) + tetenc))
    !      air%dsatdk = (610.078 * 17.27 * 237.3) / ((met%tvair-tfrz)+237.2)** 2 * &
    !                 & EXP(17.27 * (met%tvair-tfrz) / ((met%tvair-tfrz) + 237.3))  
  END SUBROUTINE define_air
END MODULE cable_air
