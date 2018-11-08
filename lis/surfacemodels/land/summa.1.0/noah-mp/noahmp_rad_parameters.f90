! ==================================================================================================
MODULE NOAHMP_RAD_PARAMETERS

    IMPLICIT NONE
 
    INTEGER I                ! loop index
    INTEGER, PARAMETER :: MSC   = 9
    INTEGER, PARAMETER :: MBAND = 2

    REAL :: ALBSAT(MSC,MBAND)   !saturated soil albedos: 1=vis, 2=nir
    REAL :: ALBDRY(MSC,MBAND)   !dry soil albedos: 1=vis, 2=nir
    REAL :: ALBICE(MBAND)       !albedo land ice: 1=vis, 2=nir
    REAL :: ALBLAK(MBAND)       !albedo frozen lakes: 1=vis, 2=nir
    REAL :: OMEGAS(MBAND)       !two-stream parameter omega for snow
    REAL :: BETADS              !two-stream parameter betad for snow
    REAL :: BETAIS              !two-stream parameter betad for snow
    REAL :: EG(2)               !emissivity

! saturated soil albedos: 1=vis, 2=nir
    DATA(ALBSAT(I,1),I=1,8)/0.15,0.11,0.10,0.09,0.08,0.07,0.06,0.05/
    DATA(ALBSAT(I,2),I=1,8)/0.30,0.22,0.20,0.18,0.16,0.14,0.12,0.10/

! dry soil albedos: 1=vis, 2=nir
    DATA(ALBDRY(I,1),I=1,8)/0.27,0.22,0.20,0.18,0.16,0.14,0.12,0.10/
    DATA(ALBDRY(I,2),I=1,8)/0.54,0.44,0.40,0.36,0.32,0.28,0.24,0.20/

! albedo land ice: 1=vis, 2=nir
    DATA (ALBICE(I),I=1,MBAND) /0.80, 0.55/

! albedo frozen lakes: 1=vis, 2=nir
    DATA (ALBLAK(I),I=1,MBAND) /0.60, 0.40/

! omega,betad,betai for snow
    DATA (OMEGAS(I),I=1,MBAND) /0.8, 0.4/
    DATA BETADS, BETAIS /0.5, 0.5/

! emissivity ground surface    
      DATA EG /0.97, 0.98/ ! 1-soil;2-lake

END MODULE NOAHMP_RAD_PARAMETERS
! ==================================================================================================


