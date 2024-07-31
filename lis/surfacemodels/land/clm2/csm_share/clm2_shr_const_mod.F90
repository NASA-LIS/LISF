!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

MODULE clm2_shr_const_mod

   use clm2_shr_kind_mod, only: CLM2_SHR_KIND_R8

   !----------------------------------------------------------------------------
   ! physical constants (all data public)
   !----------------------------------------------------------------------------
   public
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_PI     = 3.14159265358979323846_CLM2_SHR_KIND_R8  ! pi
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_CDAY   = 86400.0_CLM2_SHR_KIND_R8      ! sec in calendar day ~ sec
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_SDAY   = 86164.0_CLM2_SHR_KIND_R8      ! sec in siderial day ~ sec
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_OMEGA  = 2.0_CLM2_SHR_KIND_R8*SHR_CONST_PI/SHR_CONST_SDAY ! earth rot ~ rad sec-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_REARTH = 6.37122e6_CLM2_SHR_KIND_R8    ! radius of earth ~ m
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_G      = 9.80616_CLM2_SHR_KIND_R8      ! acceleration of gravity ~ m s-2
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_PSTD   = 101325.0_CLM2_SHR_KIND_R8     ! standard pressure ~ pascals

   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_STEBOL = 5.67e-8_CLM2_SHR_KIND_R8      ! Stefan-Boltzmann constant ~ W m-2 K-4
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_BOLTZ  = 1.38065e-23_CLM2_SHR_KIND_R8  ! Boltzmann's constant ~ J K-1 molecule-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_AVOGAD = 6.02214e26_CLM2_SHR_KIND_R8   ! Avogadro's number ~ molecules kmole-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J K-1 kmole-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_MWDAIR = 28.966_CLM2_SHR_KIND_R8       ! molecular weight dry air ~ kg kmole-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_MWWV   = 18.016_CLM2_SHR_KIND_R8       ! molecular weight water vapor
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J K-1 kg-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J K-1 kg-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_ZVIR   = (SHR_CONST_RWV/SHR_CONST_RDAIR)-1.0_CLM2_SHR_KIND_R8   ! RWV/RDAIR - 1.0
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_KARMAN = 0.4_CLM2_SHR_KIND_R8          ! Von Karman constant
 
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_TKFRZ  = 273.16_CLM2_SHR_KIND_R8       ! freezing T of fresh water ~ K (intentionally made == to TKTRIP)
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_TKTRIP = 273.16_CLM2_SHR_KIND_R8       ! triple point of fresh water ~ K

   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_RHODAIR=SHR_CONST_PSTD/ &
     (SHR_CONST_RDAIR*SHR_CONST_TKFRZ)         ! density of dry air at STP   ~ kg/m^3
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_RHOFW  = 1.000e3_CLM2_SHR_KIND_R8      ! density of fresh water ~ kg m-3
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_RHOSW  = 1.026e3_CLM2_SHR_KIND_R8      ! density of sea water ~ kg m-3
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_RHOICE = 0.917e3_CLM2_SHR_KIND_R8      ! density of ice   ~ kg m-3
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_CPDAIR = 1.00464e3_CLM2_SHR_KIND_R8    ! specific heat of dry air ~ J kg-1 K-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_CPFW   = 4.188e3_CLM2_SHR_KIND_R8      ! specific heat of fresh h2o ~ J kg-1 K-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_CPSW   = 3.996e3_CLM2_SHR_KIND_R8      ! specific heat of sea h2o ~ J kg-1 K-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_CPWV   = 1.810e3_CLM2_SHR_KIND_R8      ! specific heat of water vap ~ J kg-1 K-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_CPICE  = 2.11727e3_CLM2_SHR_KIND_R8    ! specific heat of fresh ice ~ J kg-1 K-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_LATICE = 3.337e5_CLM2_SHR_KIND_R8      ! latent heat of fusion ~ J kg-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_LATVAP = 2.501e6_CLM2_SHR_KIND_R8      ! latent heat of evaporation ~ J kg-1
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_LATSUB = SHR_CONST_LATICE + SHR_CONST_LATVAP ! latent heat of sublimation ~ J kg-1

   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_OCN_REF_SAL = 34.7_CLM2_SHR_KIND_R8    ! ocn ref salinity (psu)
   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_ICE_REF_SAL =  4.0_CLM2_SHR_KIND_R8    ! ice ref salinity (psu)

! temporary until things are converted to new names/cleaned up
!   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_TKREF = 273.16_CLM2_SHR_KIND_R8       ! conversion of K to C
!   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_ZZSICE = 0.0005_CLM2_SHR_KIND_R8      ! ice surface roughness
!   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_RHOA  = SHR_CONST_RHODAIR
!   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_RHOW  = SHR_CONST_RHOFW
!   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_RHOI  = SHR_CONST_RHOICE
!   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_CPW   = SHR_CONST_CPFW
!   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_CPAIR = SHR_CONST_CPDAIR
!   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_CPA   = SHR_CONST_CPDAIR
!   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_CW    = 0.6524_CLM2_SHR_KIND_R8
!   real(CLM2_SHR_KIND_R8),parameter :: SHR_CONST_TFRZ  = SHR_CONST_TKFRZ

END MODULE clm2_shr_const_mod
