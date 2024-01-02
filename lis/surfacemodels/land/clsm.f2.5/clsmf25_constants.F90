!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module clsmf25_constants

  ! reichle+koster, 2007
  ! reichle, 10 Oct 2008 - added "echo_catch_constants()"
  ! reichle, 28 Oct 2010 - moved DZ, SHR, PHI, FSN from subroutines gndtp0() and gndtmp()
  !                      - moved FWETL, FWETC from subroutine interc()
  !                      - renamed N_gndtmp -> N_gt
  ! reichle, 23 Nov 2010 - replaced PHIGT with POROS(N), ALHMGT with ALHM
  !                      - replaced constants with values from MAPL_Constants.F90
  !                        where possible
  ! reichle, 30 Nov 2010 - zero-diff revisions and clean-up for off-line (land-only) MERRA
  !                        replay capability
  !                         - restored PHIGT, ALHMGT 
  !                         - moved MIN_SNOW_MASS->MINSWE, DZ1MAX, and SATCAPFR to 
  !                            catch_constants
  !                         - moved "small" back from catch_constants() into snowrt()

  use clsmf25_MAPL_constants
  
  implicit none
  private
  
  public :: echo_catch_constants

  ! ---------------------------------------------------------------------------
  !
  ! physics constants

  REAL,    PARAMETER, PUBLIC :: ZERO     = 0.
  REAL,    PARAMETER, PUBLIC :: ONE      = 1.
  REAL,    PARAMETER, PUBLIC :: PIE      = MAPL_PI     ! -
  REAL,    PARAMETER, PUBLIC :: ALHE     = MAPL_ALHL   ! J/kg  @15C
  REAL,    PARAMETER, PUBLIC :: ALHM     = MAPL_ALHF   ! J/kg 
  REAL,    PARAMETER, PUBLIC :: ALHS     = MAPL_ALHS   ! J/kg
  REAL,    PARAMETER, PUBLIC :: TF       = MAPL_TICE   ! K
  REAL,    PARAMETER, PUBLIC :: STEFAN   = MAPL_STFBOL ! W/(m^2 K^4)
  REAL,    PARAMETER, PUBLIC :: RGAS     = MAPL_RGAS   ! J/(kg K)
  REAL,    PARAMETER, PUBLIC :: SHW      = MAPL_CAPWTR ! J/kg/K  spec heat of water
  REAL,    PARAMETER, PUBLIC :: SHI      = MAPL_CAPICE ! J/kg/K  spec heat of ice
  REAL,    PARAMETER, PUBLIC :: SHR      = 2400.       ! J/kg/K  spec heat of rock
  REAL,    PARAMETER, PUBLIC :: RHOW     = MAPL_RHOWTR ! kg/m^3
  REAL,    PARAMETER, PUBLIC :: GRAV     = MAPL_GRAV   ! m^2/s
  REAL,    PARAMETER, PUBLIC :: EPSILON  = MAPL_H2OMW/MAPL_AIRMW  ! -
  REAL,    PARAMETER, PUBLIC :: CPAIR    = MAPL_CP     ! J/(kg K)

  ! ---------------------------------------------------------------------------
  !
  ! constants for main Catchment model routine (catchment())

  INTEGER, PARAMETER, PUBLIC :: NTYPS    = 10     ! # vegetation types
  INTEGER, PARAMETER, PUBLIC :: N_snow   = 3      ! # layers in snow model
  INTEGER, PARAMETER, PUBLIC :: N_gt     = 6      ! # layers in ground temperature model
  INTEGER, PARAMETER, PUBLIC :: N_sm     = 3      ! # hydrological regimes considered
  
  REAL,    PARAMETER, PUBLIC :: SCONST   = 1.9E6/920.
  REAL,    PARAMETER, PUBLIC :: CSOIL_1  = 70000.
  REAL,    PARAMETER, PUBLIC :: CSOIL_2  = 200.
  
  ! ---------------------------------------------------------------------------
  !
  ! constants for snow routine (snowrt())
           
  REAL,    PARAMETER, PUBLIC :: RHOFS    = 150.  ! kg/m^3  density of fresh snow
  REAL,    PARAMETER, PUBLIC :: RHOMA    = 500.  ! kg/m^3  maximum snow density
  REAL,    PARAMETER, PUBLIC :: WEMIN    = 26.   ! kg/m^2  minimum SWE in areal fraction
  REAL,    PARAMETER, PUBLIC :: MINSWE   = 0.013 ! kg/m^2  min SWE to avoid immediate melt
  REAL,    PARAMETER, PUBLIC :: DZ1MAX   = 0.08  ! m

  ! ---------------------------------------------------------------------------
  !
  ! constants for interception routine (interc())
  
  ! Areal fraction of canopy leaves onto which precipitation falls:

  REAL,    PARAMETER, PUBLIC :: FWETL    = 0.02   ! for large-scale precipitation
  REAL,    PARAMETER, PUBLIC :: FWETC    = 0.02   ! for convective precipitation
  REAL,    PARAMETER, PUBLIC :: SATCAPFR = 0.2    ! SATCAP = SATCAPFR * LAI
  
  ! ---------------------------------------------------------------------------
  !
  ! constants for ground temperature routine (gndtp0() and gndtmp())

  REAL,    PARAMETER, DIMENSION(N_gt), PUBLIC :: DZGT = &  ! m  layer depths 
       (/ 0.0988, 0.1952, 0.3859, 0.7626, 1.5071, 10.0 /)

  ! PHIGT and ALHMGT are needed for backward compatibility with 
  !  off-line (land-only) MERRA replay:
  ! 
  ! PHIGT = porosity used in gndtp0() and gndtmp() 
  !         if neg,  POROS(n) from soil moisture submodel will be used
  ! 
  !               |   PHIGT      ALHMGT
  ! ------------------------------------------------
  !  MERRA        |      0.45    3.34e+5   
  !  Fortuna-2_3  |  -9999.       ALHM

  REAL,    PARAMETER, PUBLIC :: PHIGT   = -9999.     
  REAL,    PARAMETER, PUBLIC :: ALHMGT  = ALHM

  
contains
  
  subroutine echo_catch_constants(logunit)
    
    ! reichle, 10 Oct 2008
    
    implicit none
    
    integer, intent(in) :: logunit
    
    write (logunit,*)
    write (logunit,*) '-----------------------------------------------------------'
    write (logunit,*)
    write (logunit,*) 'echo_catch_constants():'
    write (logunit,*)
    write (logunit,*) 'ZERO     = ', ZERO     
    write (logunit,*) 'ONE      = ', ONE      
    write (logunit,*) 'PIE      = ', PIE      
    write (logunit,*) 'ALHE     = ', ALHE     
    write (logunit,*) 'ALHM     = ', ALHM     
    write (logunit,*) 'ALHS     = ', ALHS     
    write (logunit,*) 'TF       = ', TF    
    write (logunit,*) 'STEFAN   = ', STEFAN
    write (logunit,*) 'RGAS     = ', RGAS     
    write (logunit,*) 'SHW      = ', SHW      
    write (logunit,*) 'SHI      = ', SHI      
    write (logunit,*) 'SHR      = ', SHR    
    write (logunit,*) 'RHOW     = ', RHOW     
    write (logunit,*) 'GRAV     = ', GRAV     
    write (logunit,*) 'EPSILON  = ', EPSILON  
    write (logunit,*) 'CPAIR    = ', CPAIR    
    write (logunit,*)
    write (logunit,*) 'NTYPS    = ', NTYPS     
    write (logunit,*) 'N_snow   = ', N_snow   
    write (logunit,*) 'N_gt     = ', N_gt 
    write (logunit,*) 'N_sm     = ', N_sm     
    write (logunit,*)
    write (logunit,*) 'SCONST   = ', SCONST   
    write (logunit,*) 'CSOIL_1  = ', CSOIL_1    
    write (logunit,*) 'CSOIL_2  = ', CSOIL_2    
    write (logunit,*)
    write (logunit,*) 'RHOFS    = ', RHOFS    
    write (logunit,*) 'RHOMA    = ', RHOMA    
    write (logunit,*) 'WEMIN    = ', WEMIN    
    write (logunit,*) 'MINSWE   = ', MINSWE
    write (logunit,*) 'DZ1MAX   = ', DZ1MAX   
    write (logunit,*) 
    write (logunit,*) 'FWETL    = ', FWETL     
    write (logunit,*) 'FWETC    = ', FWETC     
    write (logunit,*) 'SATCAPFR = ', SATCAPFR
    write (logunit,*) 
    write (logunit,*) 'DZGT     = ', DZGT 
    write (logunit,*) 'PHIGT    = ', PHIGT
    write (logunit,*) 'ALHMGT   = ', ALHMGT 
    write (logunit,*) 
    write (logunit,*) 'end echo_catch_constants()'
    write (logunit,*)
    write (logunit,*) '-----------------------------------------------------------'
    write (logunit,*)
    
  end subroutine echo_catch_constants
  
end module clsmf25_constants

! Previous LIS-Catchment values:
!  REAL, PARAMETER :: ALHE = 2.4548E6
!  REAL, PARAMETER :: ALHS = 2.8368E6
!  REAL, PARAMETER :: STEFAN = 5.669E-8
!  REAL, PARAMETER :: RGAS = .286*1003.5
!  REAL, PARAMETER :: SHW = 4185.
!  REAL, PARAMETER :: SHI = 2060.
!  REAL, PARAMETER :: GRAV = 9.81
!  REAL, PARAMETER :: CPAIR = 1010.
!  REAL, PARAMETER :: RHOMA  = 700.   !  [KG/M^3] maximum snow density
!  REAL, PARAMETER :: WEMIN  = 13.    !  [KG/M^2] minimum SWE in areal fraction
!  REAL, PARAMETER :: DZ1MAX = 0.05    ! m

