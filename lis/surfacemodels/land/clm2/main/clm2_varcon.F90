!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

module clm2_varcon

!----------------------------------------------------------------------- 
! 
! Purpose: 
! module for land model constants 

! Method: 
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: clm2_varcon.F90,v 1.7 2004/11/24 22:57:10 jim Exp $
!-----------------------------------------------------------------------

!<debug>
!  use LIS_precisionMod
  use LIS_precisionMod, only : r8
!<debug>
  use clm2_shr_const_mod, only: SHR_CONST_G,SHR_CONST_STEBOL,SHR_CONST_KARMAN,     &
                           SHR_CONST_RWV,SHR_CONST_RDAIR,SHR_CONST_CPFW,      &
                           SHR_CONST_CPICE,SHR_CONST_CPDAIR,SHR_CONST_LATVAP, &
                           SHR_CONST_LATSUB,SHR_CONST_LATICE,SHR_CONST_RHOFW, &
                           SHR_CONST_RHOICE,SHR_CONST_TKFRZ,SHR_CONST_REARTH
  use clm2_varpar, only : numcol, numrad
  implicit none

  PRIVATE
!------------------------------------------------------------------
! Initialize physical constants
!------------------------------------------------------------------

  real(r8), public  :: grav   = SHR_CONST_G      !gravity constant [m s-12]
  real(r8), public  :: sb     = SHR_CONST_STEBOL !stefan-boltzmann constant  [W m-2/K4]
  real(r8), public  :: vkc    = SHR_CONST_KARMAN !von Karman constant [-]
  real(r8), public  :: rwat   = SHR_CONST_RWV    !gas constant for water vapor [J/(kg K)]
  real(r8), public  :: rair   = SHR_CONST_RDAIR  !gas constant for dry air [J kg-1/K]
  real(r8), public  :: roverg = SHR_CONST_RWV/SHR_CONST_G*1000. !Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  real(r8), public  :: cpliq  = SHR_CONST_CPFW   !Specific heat of water [J kg-1-K]
  real(r8), public  :: cpice  = SHR_CONST_CPICE  !Specific heat of ice [J kg-1-K]
  real(r8), public  :: cpair  = SHR_CONST_CPDAIR !specific heat of dry air [J kg-1/K]
  real(r8), public  :: hvap   = SHR_CONST_LATVAP !Latent heat of evap for water [J kg-1]
  real(r8), public  :: hsub   = SHR_CONST_LATSUB !Latent heat of sublimation    [J kg-1]
  real(r8), public  :: hfus   = SHR_CONST_LATICE !Latent heat of fusion for ice [J kg-1]
  real(r8), public  :: denh2o = SHR_CONST_RHOFW  !density of liquid water [kg/m3]
  real(r8), public  :: denice = SHR_CONST_RHOICE !density of ice [kg/m3]
  real(r8), public  :: tkair  = 0.023     !thermal conductivity of air   [W/m/k]
  real(r8), public  :: tkice  = 2.290     !thermal conductivity of ice   [W/m/k]
  real(r8), public  :: tkwat  = 0.6       !thermal conductivity of water [W/m/k]
  real(r8), public  :: tfrz   = SHR_CONST_TKFRZ  !freezing temperature [K]
  real(r8), public  :: tcrit  = 2.5       !critical temperature to determine rain or snow
  real(r8), public  :: po2    = 0.209     !constant atmospheric partial pressure  O2 (mol/mol)
  real(r8), public  :: pco2   = 355.e-06  !constant atmospheric partial pressure CO2 (mol/mol)

  real(r8), public  :: bdsno = 250.       !bulk density snow (kg/m**3)

  real(r8), public  :: re = SHR_CONST_REARTH*0.001 !radius of earth (km)

  real(r8), public, parameter :: spval = 1.e36  !special value for missing data (ocean)

!------------------------------------------------------------------
! Initialize water type constants
!------------------------------------------------------------------

! "water" types 
!   1     soil
!   2     land ice (glacier)
!   3     deep lake
!   4     shallow lake
!   5     wetland: swamp, marsh, etc

  integer, public :: istsoil = 1  !soil         "water" type
  integer, public :: istice  = 2  !land ice     "water" type
  integer, public :: istdlak = 3  !deep lake    "water" type
  integer, public :: istslak = 4  !shallow lake "water" type
  integer, public :: istwet  = 5  !wetland      "water" type

!------------------------------------------------------------------
! Initialize miscellaneous radiation constants
!------------------------------------------------------------------

  integer, private :: i  ! loop index

! saturated soil albedos for 8 color classes: 1=vis, 2=nir

  real(r8), public :: albsat(numcol,numrad) !wet soil albedo by color class and waveband
  data(albsat(i,1),i=1,8)/0.12,0.11,0.10,0.09,0.08,0.07,0.06,0.05/
  data(albsat(i,2),i=1,8)/0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10/

! dry soil albedos for 8 color classes: 1=vis, 2=nir 

  real(r8), public :: albdry(numcol,numrad) !dry soil albedo by color class and waveband
  data(albdry(i,1),i=1,8)/0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10/
  data(albdry(i,2),i=1,8)/0.48,0.44,0.40,0.36,0.32,0.28,0.24,0.20/

! albedo land ice: 1=vis, 2=nir

  real(r8), public :: albice(numrad)        !albedo land ice by waveband
  data (albice(i),i=1,numrad) /0.80, 0.55/

! albedo frozen lakes: 1=vis, 2=nir 

  real(r8), public :: alblak(numrad)        !albedo frozen lakes by waveband
  data (alblak(i),i=1,numrad) /0.60, 0.40/

! omega,betad,betai for snow 

  real(r8), public :: betads  = 0.5       !two-stream parameter betad for snow
  real(r8), public :: betais  = 0.5       !two-stream parameter betai for snow
  real(r8), public :: omegas(numrad)      !two-stream parameter omega for snow by band
  data (omegas(i),i=1,numrad) /0.8, 0.4/
  integer, public  :: iyear_AD ! to simulate above earth's orbital parameters for
  logical, public :: doalb
  real(r8), public :: eccen   ! Earth's eccentricity factor
  real(r8), public :: obliq   ! Earth's obliquity angle
  real(r8), public :: mvelpp   ! Earth's moving vernal equinoz at perhelion
!===  Orbital information after call to routine clm2_shr_orbit_params
  real(r8), public :: obliqr  ! Earth's obliquity in radians
  real(r8), public :: lambm0  ! Mean longitude (radians) of perihelion at the vernal equinox
  logical :: log_print
end module clm2_varcon
