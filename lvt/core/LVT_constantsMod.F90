!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LVT_constantsMod
!BOP
!
!  !MODULE: LVT_constantsMod
! 
!  !DESCRIPTION: 
!   The code in this file provides values of physical constants for
!   consistent use across different components. 
!   
!  !REVISION HISTORY: 
!  17 Feb 2004    Sujay Kumar  Initial Specification
!
!EOP
!BOC
!----------------------------------------------------------------------------
! physical constants (all data public)
!----------------------------------------------------------------------------
   public
!   real,parameter :: CONST_PI     = 3.14159265358979323846  ! pi
   real,parameter :: LVT_CONST_PI     = 3.14159265   ! pi
   real,parameter :: LVT_CONST_CDAY   = 86400.0      ! sec in calendar day ~ sec
   real,parameter :: LVT_CONST_SDAY   = 86164.0      ! sec in siderial day ~ sec
   real,parameter :: LVT_CONST_OMEGA  = 2.0*LVT_CONST_PI/LVT_CONST_SDAY ! earth rot ~ rad/sec
   real,parameter :: LVT_CONST_REARTH = 6.37122e6    ! radius of earth ~ m
   real,parameter :: LVT_CONST_G      = 9.80616      ! acceleration of gravity ~ m/s^2
   real,parameter :: LVT_CONST_PSTD   = 101325.0     ! standard pressure ~ pascals

   real,parameter :: LVT_CONST_STEBOL = 5.67e-8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
   real,parameter :: LVT_CONST_BOLTZ  = 1.38065e-23  ! Boltzmann's constant ~ J/K/molecule
   real,parameter :: LVT_CONST_AVOGAD = 6.02214e26   ! Avogadro's number ~ molecules/kmole
   real,parameter :: LVT_CONST_RGAS   = LVT_CONST_AVOGAD*LVT_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
   real,parameter :: LVT_CONST_MWDAIR = 28.966       ! molecular weight dry air ~ kg/kmole
   real,parameter :: LVT_CONST_MWWV   = 18.016       ! molecular weight water vapor
   real,parameter :: LVT_CONST_RDAIR  = LVT_CONST_RGAS/LVT_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
   real,parameter :: LVT_CONST_RWV    = LVT_CONST_RGAS/LVT_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg
   real,parameter :: LVT_CONST_ZVIR   = (LVT_CONST_RWV/LVT_CONST_RDAIR)-1.0   ! RWV/RDAIR - 1.0
   real,parameter :: LVT_CONST_KARMAN = 0.4          ! Von Karman constant
 
   real,parameter :: LVT_CONST_TKFRZ  = 273.16       ! freezing T of fresh water ~ K (intentionally made == to TKTRIP)
   real,parameter :: LVT_CONST_TKTRIP = 273.16       ! triple point of fresh water ~ K

   real,parameter :: LVT_CONST_RHODAIR=LVT_CONST_PSTD/ &
     (LVT_CONST_RDAIR*LVT_CONST_TKFRZ)         ! density of dry air at STP   ~ kg/m^3
   real,parameter :: LVT_CONST_RHOFW  = 1.000e3      ! density of fresh water ~ kg/m^3
   real,parameter :: LVT_CONST_RHOSW  = 1.026e3      ! density of sea water ~ kg/m^3
   real,parameter :: LVT_CONST_RHOICE = 0.917e3      ! density of ice   ~ kg/m^3
   real,parameter :: LVT_CONST_CPDAIR = 1.00464e3    ! specific heat of dry air ~ J/kg/K
   real,parameter :: LVT_CONST_CPFW   = 4.188e3      ! specific heat of fresh h2o ~ J/kg/K
   real,parameter :: LVT_CONST_CPSW   = 3.996e3      ! specific heat of sea h2o ~ J/kg/K
   real,parameter :: LVT_CONST_CPWV   = 1.810e3      ! specific heat of water vap ~ J/kg/K
   real,parameter :: LVT_CONST_CPICE  = 2.11727e3    ! specific heat of fresh ice ~ J/kg/K
   real,parameter :: LVT_CONST_LATICE = 3.337e5      ! latent heat of fusion ~ J/kg
   real,parameter :: LVT_CONST_LATVAP = 2.501e6      ! latent heat of evaporation ~ J/kg
   real,parameter :: LVT_CONST_LATSUB = LVT_CONST_LATICE + LVT_CONST_LATVAP ! latent heat of sublimation ~ J/kg

   real,parameter :: LVT_CONST_OCN_REF_SAL = 34.7    ! ocn ref salinity (psu)
   real,parameter :: LVT_CONST_ICE_REF_SAL =  4.0    ! ice ref salinity (psu)
   real,parameter :: LVT_CONST_SOLAR  =  1369.2    !Solar constant (W/m2)

   real, parameter :: LVT_MS2KMDAY = 86.4 ! m/s to km/day
   real, parameter :: LVT_CONST_ANGULAR_VELOCITY = 0.2618 ! Angular velocity of earth's rotation (radian/hr)

!EOC
 end module LVT_constantsMod
