MODULE cable_sli_numbers

  USE cable_dimensions,  ONLY: r_2, i_d
  USE cable_physical_constants, ONLY: tfrz, grav, rmh2o, rhow, rgas, capp, tetena, tetenb, tetenc

  IMPLICIT NONE

  ! define types
  TYPE vars_met
     REAL(r_2) :: Ta, rha, rbw, rbh, rrc, ra, rs, Rn, u, Da, cva, civa, ha, phiva, qevappot
  END TYPE vars_met

  TYPE vars
     INTEGER(i_d) :: isat
     REAL(r_2)    :: h, phi, phiS, K, KS, Dv, cvsat, rh, phiv, phivS, kH
     REAL(r_2)    :: kE, kth, csoil, eta_th, hS, rhS, sl, cv, cvsatT, cvS, kv
     INTEGER(i_d) :: iice
     REAL(r_2)    :: thetai, thetal, phiT, KT, lambdav, lambdaf
     REAL(r_2)    :: he, phie, Ksat ! air-entry potential values (different to core sp params for frozen soil)
     REAL(r_2)    :: dthetaldT
     REAL(r_2)    :: Tfrz ! freezing point (deg C) depends on csol and S
     REAL(r_2)    :: csoileff
  END TYPE vars

  TYPE vars_aquifer
     INTEGER(i_d) :: isat
     REAL(r_2)    :: zdelta, zsoil, zzero, K, Wa, discharge, f, Rsmax, Sy
  END TYPE vars_aquifer

  TYPE params
     REAL(r_2) :: the, thre, he, lam, Ke, eta,thr
     REAL(r_2) :: KSe, phie, phiSe, rho, thw, kd, css, clay, tortuosity
     REAL(r_2) :: isahorizon, isbhorizon
  END TYPE params

  TYPE rapointer
     REAL(r_2), DIMENSION(:), POINTER :: p
  END TYPE rapointer

  TYPE solve_type ! for energy and moisture balance in rh0_sol, etc.
     INTEGER(i_d) :: k
     REAL(r_2)    :: T1, Ta, cva, Rnet, hr1, hra, Dv, gv, gh, Dh, dz, phie, he, K1, eta,lambda, Ks, lambdav
  END TYPE solve_type

  ! define some numbers
  REAL(r_2), PARAMETER :: zero      = 0.0
  REAL(r_2), PARAMETER :: half      = 0.5
  REAL(r_2), PARAMETER :: one       = 1.0
  REAL(r_2), PARAMETER :: two       = 2.0
  REAL(r_2), PARAMETER :: four      = 4.0
  REAL(r_2), PARAMETER :: thousand  = 1000.
  REAL(r_2), PARAMETER :: e1        = 1.e-1
  REAL(r_2), PARAMETER :: e2        = 1.e-2
  REAL(r_2), PARAMETER :: e3        = 1.e-3
  REAL(r_2), PARAMETER :: e4        = 1.e-4
  REAL(r_2), PARAMETER :: e5        = 1.e-5
  REAL(r_2), PARAMETER :: e6        = 1.e-6
  REAL(r_2), PARAMETER :: e7        = 1.e-7

  ! define some constants
  REAL(r_2), PARAMETER :: Tzero     = tfrz      ! Celcius -> Kelvin
  REAL(r_2), PARAMETER :: gravity   = grav      ! gravitation constant [m/s2]
  REAL(r_2), PARAMETER :: Mw        = rmh2o     ! weight of 1 mol of water [kg]
  REAL(r_2), PARAMETER :: cpa       = capp      ! specific heat capacity of dry air at 0-40 degC [J/kgK]
  REAL(r_2), PARAMETER :: esata     = tetena*100.0_r_2 ! constants for saturated vapour pressure calculation
  REAL(r_2), PARAMETER :: esatb     = tetenb    ! %
  REAL(r_2), PARAMETER :: esatc     = tetenc    ! %

  REAL(r_2), PARAMETER :: rlambda   = 2.442e6   ! latent heat of condensation at 25 degC [J/kg]
  REAL(r_2), PARAMETER :: lambdaf   = 335000.   ! latent heat of fusion (J kg-1)
  REAL(r_2), PARAMETER :: lambdas   = 2835000.  ! latent heat of sublimation (J kg-1)
  REAL(r_2), PARAMETER :: Dva       = 2.17e-5   ! vapour diffusivity of water in air at 0 degC [m2/s]
  !REAL(r_2), PARAMETER :: rhow      = 1000.    ! denisty of water [kg/m3]
  REAL(r_2), PARAMETER :: rhoi      = 920.      ! density of ice (kg m-3) 

  REAL(r_2), PARAMETER :: rhoa      = 1.184     ! denisty of dry air at std (25 degC) [kg/m3]
  REAL(r_2), PARAMETER :: rhocp     = 1189.8    ! cpa*rhoa at std (25 degC) [J/m3K]

  REAL(r_2), PARAMETER :: esata_ice = 611.2   ! constants for saturated vapour pressure calculation over ice (WMO, 2008)
  REAL(r_2), PARAMETER :: esatb_ice = 22.46   ! %
  REAL(r_2), PARAMETER :: esatc_ice = 272.62  ! %
  REAL(r_2), PARAMETER :: csice     = 2.100e3 ! specific heat capacity for ice
  !REAL(r_2), PARAMETER :: csice     = 4.218e3 ! specific heat capacity for ice
  REAL(r_2), PARAMETER :: cswat     = 4.218e3 ! specific heat capacity for water

!!!MC!!!
  REAL(r_2), PARAMETER :: kw        = 0.58    ! dito

  ! numerical limits
  REAL(r_2), PARAMETER :: dSfac     = 1.25
  REAL(r_2), PARAMETER :: dpmaxr    = 0.5
  REAL(r_2), PARAMETER :: h0min     = -2.e-4
  REAL(r_2), PARAMETER :: Smax      = 1.001
  REAL(r_2), PARAMETER :: dh0max    = 0.0001
  REAL(r_2), PARAMETER :: SLmax     = 1.01
  REAL(r_2), PARAMETER :: SLmin     = 0.001
  REAL(r_2), PARAMETER :: h0max     = 0.01e-0
  REAL(r_2), PARAMETER :: qprecmax  = 1.0e10
  REAL(r_2), PARAMETER :: dSmax     = 0.5
  REAL(r_2), PARAMETER :: dSmaxr    = 0.1
  REAL(r_2), PARAMETER :: dtmax     = 3600.
  REAL(r_2), PARAMETER :: dtmin     = 0.01
  REAL(r_2), PARAMETER :: dsmmax    = 1.0
  REAL(r_2), PARAMETER :: dTsoilmax = 5.0
  REAL(r_2), PARAMETER :: dTLmax    = 0.5
  REAL(r_2), PARAMETER :: tol_dthetaldT = 1.e-6_r_2
  INTEGER(i_d), PARAMETER ::nsteps_ice_max = 20
  !MC-ToDo! Identify why smaller time-steps are needed for isotopes
  ! With isotopes  REAL(r_2), PARAMETER :: dSmax=0.1_r_2, dSmaxr=0.1_r_2, dtmax=0.05_r_2*3600., 
  !                                        dtmin=0.0_r_2, dsmmax=1.0_r_2
  ! With isotopes   REAL(r_2), PARAMETER ::  dTsoilmax=1.0_r_2, dTLmax=1.0_r_2
  REAL(r_2), PARAMETER :: gf        = 1.0    ! gravity factor for flow direction (usually 1 for straight down).
  REAL(r_2), PARAMETER :: hmin      = -1.0e6 ! minimum matric head h (used by MF).
  REAL(r_2), PARAMETER :: csol      = 0.0    ! solute concentration (mol kg-1 soil water)
  REAL(r_2), PARAMETER :: rhmin     = 0.05   ! minimum relative humidity in soil and litter

  ! boundary conditions
  REAL(r_2),         PARAMETER :: hbot  = 0.0
  CHARACTER(LEN=20), PARAMETER :: botbc = "free drainage"
  !CHARACTER(LEN=20), PARAMETER :: botbc = "zero flux"
  !CHARACTER(LEN=20), PARAMETER :: botbc = "aquifer"
  !CHARACTER(LEN=20), PARAMETER :: botbc = "constant head"
  !CHARACTER(LEN=20), PARAMETER :: botbc = "seepage"

END MODULE cable_sli_numbers
