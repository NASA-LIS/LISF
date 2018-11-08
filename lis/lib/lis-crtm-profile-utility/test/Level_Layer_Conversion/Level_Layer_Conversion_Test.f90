!
! Level_Layer_Conversion_Test
!
! Program to test the routines in the Level_Layer_Conversion module
! of the Profile_Utility library.
!
! Note that the Create_Sublevels and Integrate_Sublevels routines are
! based on and adapted from the UMBC KLAYERS code.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 01-Dec-2000
!                       paul.vandelst@ssec.wisc.edu
!

PROGRAM Level_Layer_Conversion_Test

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module usage
  USE Profile_Utility
  ! Disable all implicit typing
  IMPLICIT NONE

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'Level_Layer_Conversion_Test'
  CHARACTER(*), PARAMETER :: PROGRAM_RCS_ID = &
    '$Id: Level_Layer_Conversion_Test.f90 1958 2008-03-18 16:35:36Z paul.vandelst@noaa.gov $'
  ! Output filename
  CHARACTER(*), PARAMETER :: FILENAME = 'Profile_Conversion.asc'
  ! The dimensions
  INTEGER, PARAMETER :: N_LEVELS    = 50
  INTEGER, PARAMETER :: N_LAYERS    = N_LEVELS - 1
  INTEGER, PARAMETER :: N_ABSORBERS = 6
  INTEGER, PARAMETER :: N_PER_LAYER = 10
  INTEGER, PARAMETER :: N_SUBLEVELS = ( N_LAYERS * N_PER_LAYER ) + 1
  ! Location of water vapour in absorber array
  INTEGER, PARAMETER :: H2O_J_INDEX = 1

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  INTEGER  :: FileID
  INTEGER  :: j, k
  INTEGER  :: IO_Status
  INTEGER  :: Error_Status
  ! Profile data
  REAL(fp) :: Pressure(N_LEVELS)
  REAL(fp) :: Temperature(N_LEVELS)
  REAL(fp) :: Absorber(N_LEVELS,N_ABSORBERS)
  ! Output from PPMV_to_PP
  REAL(fp) :: H2O_pp(N_LEVELS)
  REAL(fp) :: Sublevel_H2O_pp(N_SUBLEVELS)
  ! Outputs from Geopotential_Height
  REAL(fp) :: z(N_LEVELS)
  REAL(fp) :: Sublevel_z(N_SUBLEVELS)
  ! Outputs from Effective_Layer_TP
  REAL(fp) :: Effective_Layer_p(N_LAYERS)
  REAL(fp) :: Effective_Layer_t(N_LAYERS)
  ! Outputs from Create_Sublevels
  REAL(fp) :: Sublevel_p(N_SUBLEVELS) 
  REAL(fp) :: Sublevel_t(N_SUBLEVELS) 
  REAL(fp) :: Sublevel_aa(N_SUBLEVELS,N_ABSORBERS)
  ! Outputs from Integrate_Sublevels
  REAL(fp) :: Layer_p(N_LAYERS)            
  REAL(fp) :: Layer_t(N_LAYERS)            
  REAL(fp) :: Layer_aa(N_LAYERS,N_ABSORBERS)

  ! Output program header
  ! ---------------------
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to test the Profile_Utility Level_Layer_Conversion module '//&
                        'routines. Note that the Create_Sublevels and Integrate_Sublevels '//&
                        'routines are based on and adapted from the UMBC KLAYERS code.', &
                        '$Revision: 1958 $' )

  ! Load the US Std Atm profile data
  ! --------------------------------
  Pressure = (/ &
    1.013e+03_fp, 8.988e+02_fp, 7.950e+02_fp, 7.012e+02_fp, 6.166e+02_fp, &
    5.405e+02_fp, 4.722e+02_fp, 4.111e+02_fp, 3.565e+02_fp, 3.080e+02_fp, &
    2.650e+02_fp, 2.270e+02_fp, 1.940e+02_fp, 1.658e+02_fp, 1.417e+02_fp, &
    1.211e+02_fp, 1.035e+02_fp, 8.850e+01_fp, 7.565e+01_fp, 6.467e+01_fp, &
    5.529e+01_fp, 4.729e+01_fp, 4.047e+01_fp, 3.467e+01_fp, 2.972e+01_fp, &
    2.549e+01_fp, 1.743e+01_fp, 1.197e+01_fp, 8.010e+00_fp, 5.746e+00_fp, &
    4.150e+00_fp, 2.871e+00_fp, 2.060e+00_fp, 1.491e+00_fp, 1.090e+00_fp, &
    7.978e-01_fp, 4.250e-01_fp, 2.190e-01_fp, 1.090e-01_fp, 5.220e-02_fp, &
    2.400e-02_fp, 1.050e-02_fp, 4.460e-03_fp, 1.840e-03_fp, 7.600e-04_fp, &
    3.200e-04_fp, 1.450e-04_fp, 7.100e-05_fp, 4.010e-05_fp, 2.540e-05_fp /)
  Temperature = (/ &
    288.20_fp, 281.70_fp, 275.20_fp, 268.70_fp, 262.20_fp, &
    255.70_fp, 249.20_fp, 242.70_fp, 236.20_fp, 229.70_fp, &
    223.30_fp, 216.80_fp, 216.70_fp, 216.70_fp, 216.70_fp, &
    216.70_fp, 216.70_fp, 216.70_fp, 216.70_fp, 216.70_fp, &
    216.70_fp, 217.60_fp, 218.60_fp, 219.60_fp, 220.60_fp, &
    221.60_fp, 224.00_fp, 226.50_fp, 230.00_fp, 236.50_fp, &
    242.90_fp, 250.40_fp, 257.30_fp, 264.20_fp, 270.60_fp, &
    270.70_fp, 260.80_fp, 247.00_fp, 233.30_fp, 219.60_fp, &
    208.40_fp, 198.60_fp, 188.90_fp, 186.90_fp, 188.40_fp, &
    195.10_fp, 208.80_fp, 240.00_fp, 300.00_fp, 360.00_fp /)
  ! Water vapour
  Absorber(:,1) = (/ &
    7.745e+03_fp, 6.071e+03_fp, 4.631e+03_fp, 3.182e+03_fp, 2.158e+03_fp, &
    1.397e+03_fp, 9.254e+02_fp, 5.720e+02_fp, 3.667e+02_fp, 1.583e+02_fp, &
    6.996e+01_fp, 3.613e+01_fp, 1.906e+01_fp, 1.085e+01_fp, 5.927e+00_fp, &
    5.000e+00_fp, 3.950e+00_fp, 3.850e+00_fp, 3.825e+00_fp, 3.850e+00_fp, &
    3.900e+00_fp, 3.975e+00_fp, 4.065e+00_fp, 4.200e+00_fp, 4.300e+00_fp, &
    4.425e+00_fp, 4.575e+00_fp, 4.725e+00_fp, 4.825e+00_fp, 4.900e+00_fp, &
    4.950e+00_fp, 5.025e+00_fp, 5.150e+00_fp, 5.225e+00_fp, 5.250e+00_fp, &
    5.225e+00_fp, 5.100e+00_fp, 4.750e+00_fp, 4.200e+00_fp, 3.500e+00_fp, &
    2.825e+00_fp, 2.050e+00_fp, 1.330e+00_fp, 8.500e-01_fp, 5.400e-01_fp, &
    4.000e-01_fp, 3.400e-01_fp, 2.800e-01_fp, 2.400e-01_fp, 2.000e-01_fp /)
  ! Ozone
  Absorber(:,2) = (/ &
    2.660e-02_fp, 2.931e-02_fp, 3.237e-02_fp, 3.318e-02_fp, 3.387e-02_fp, &
    3.768e-02_fp, 4.112e-02_fp, 5.009e-02_fp, 5.966e-02_fp, 9.168e-02_fp, &
    1.313e-01_fp, 2.149e-01_fp, 3.095e-01_fp, 3.846e-01_fp, 5.030e-01_fp, &
    6.505e-01_fp, 8.701e-01_fp, 1.187e+00_fp, 1.587e+00_fp, 2.030e+00_fp, &
    2.579e+00_fp, 3.028e+00_fp, 3.647e+00_fp, 4.168e+00_fp, 4.627e+00_fp, &
    5.118e+00_fp, 5.803e+00_fp, 6.553e+00_fp, 7.373e+00_fp, 7.837e+00_fp, &
    7.800e+00_fp, 7.300e+00_fp, 6.200e+00_fp, 5.250e+00_fp, 4.100e+00_fp, &
    3.100e+00_fp, 1.800e+00_fp, 1.100e+00_fp, 7.000e-01_fp, 3.000e-01_fp, &
    2.500e-01_fp, 3.000e-01_fp, 5.000e-01_fp, 7.000e-01_fp, 7.000e-01_fp, &
    4.000e-01_fp, 2.000e-01_fp, 5.000e-02_fp, 5.000e-03_fp, 5.000e-04_fp /)
  ! Nitrous oxide
  Absorber(:,3) = (/ &
    3.200e-01_fp, 3.200e-01_fp, 3.200e-01_fp, 3.200e-01_fp, 3.200e-01_fp, &
    3.200e-01_fp, 3.200e-01_fp, 3.200e-01_fp, 3.200e-01_fp, 3.195e-01_fp, &
    3.179e-01_fp, 3.140e-01_fp, 3.095e-01_fp, 3.048e-01_fp, 2.999e-01_fp, &
    2.944e-01_fp, 2.877e-01_fp, 2.783e-01_fp, 2.671e-01_fp, 2.527e-01_fp, &
    2.365e-01_fp, 2.194e-01_fp, 2.051e-01_fp, 1.967e-01_fp, 1.875e-01_fp, &
    1.756e-01_fp, 1.588e-01_fp, 1.416e-01_fp, 1.165e-01_fp, 9.275e-02_fp, &
    6.693e-02_fp, 4.513e-02_fp, 2.751e-02_fp, 1.591e-02_fp, 9.378e-03_fp, &
    4.752e-03_fp, 3.000e-03_fp, 2.065e-03_fp, 1.507e-03_fp, 1.149e-03_fp, &
    8.890e-04_fp, 7.056e-04_fp, 5.716e-04_fp, 4.708e-04_fp, 3.932e-04_fp, &
    3.323e-04_fp, 2.837e-04_fp, 2.443e-04_fp, 2.120e-04_fp, 1.851e-04_fp /)
  ! Carbon monoxide
  Absorber(:,4) = (/ &
    1.500e-01_fp, 1.450e-01_fp, 1.399e-01_fp, 1.349e-01_fp, 1.312e-01_fp, &
    1.303e-01_fp, 1.288e-01_fp, 1.247e-01_fp, 1.185e-01_fp, 1.094e-01_fp, &
    9.962e-02_fp, 8.964e-02_fp, 7.814e-02_fp, 6.374e-02_fp, 5.025e-02_fp, &
    3.941e-02_fp, 3.069e-02_fp, 2.489e-02_fp, 1.966e-02_fp, 1.549e-02_fp, &
    1.331e-02_fp, 1.232e-02_fp, 1.232e-02_fp, 1.307e-02_fp, 1.400e-02_fp, &
    1.498e-02_fp, 1.598e-02_fp, 1.710e-02_fp, 1.850e-02_fp, 2.009e-02_fp, &
    2.220e-02_fp, 2.497e-02_fp, 2.824e-02_fp, 3.241e-02_fp, 3.717e-02_fp, &
    4.597e-02_fp, 6.639e-02_fp, 1.073e-01_fp, 1.862e-01_fp, 3.059e-01_fp, &
    6.375e-01_fp, 1.497e+00_fp, 3.239e+00_fp, 5.843e+00_fp, 1.013e+01_fp, &
    1.692e+01_fp, 2.467e+01_fp, 3.356e+01_fp, 4.148e+01_fp, 5.000e+01_fp /)
  ! Methane
  Absorber(:,5) = (/ &
    1.700e+00_fp, 1.700e+00_fp, 1.700e+00_fp, 1.700e+00_fp, 1.700e+00_fp, &
    1.700e+00_fp, 1.700e+00_fp, 1.699e+00_fp, 1.697e+00_fp, 1.693e+00_fp, &
    1.685e+00_fp, 1.675e+00_fp, 1.662e+00_fp, 1.645e+00_fp, 1.626e+00_fp, &
    1.605e+00_fp, 1.582e+00_fp, 1.553e+00_fp, 1.521e+00_fp, 1.480e+00_fp, &
    1.424e+00_fp, 1.355e+00_fp, 1.272e+00_fp, 1.191e+00_fp, 1.118e+00_fp, &
    1.055e+00_fp, 9.870e-01_fp, 9.136e-01_fp, 8.300e-01_fp, 7.460e-01_fp, &
    6.618e-01_fp, 5.638e-01_fp, 4.614e-01_fp, 3.631e-01_fp, 2.773e-01_fp, &
    2.100e-01_fp, 1.650e-01_fp, 1.500e-01_fp, 1.500e-01_fp, 1.500e-01_fp, &
    1.500e-01_fp, 1.500e-01_fp, 1.500e-01_fp, 1.400e-01_fp, 1.300e-01_fp, &
    1.200e-01_fp, 1.100e-01_fp, 9.500e-02_fp, 6.000e-02_fp, 3.000e-02_fp /)
  ! Nitrogen
  Absorber(:,6) = (/ &
    7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, &
    7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, &
    7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, &
    7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, &
    7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, &
    7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, &
    7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, &
    7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, &
    7.81e+05_fp, 7.81e+05_fp, 7.81e+05_fp, 7.80e+05_fp, 7.79e+05_fp, &
    7.77e+05_fp, 7.74e+05_fp, 7.70e+05_fp, 7.65e+05_fp, 7.60e+05_fp /)

  ! Open the output data file
  ! -------------------------
  FileID = Get_Lun()
  IF ( FileID < 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error obtaining output file unit number', &
                          FAILURE )
    STOP
  END IF
  OPEN( FileID,FILE  =FILENAME, &
               STATUS='REPLACE', &
               FORM  ='FORMATTED', &
               ACCESS='SEQUENTIAL', &
               IOSTAT=IO_Status )
  IF ( IO_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error opening output file '//FILENAME, &
                          FAILURE )
    STOP
  END IF


  ! ======================================================
  ! CREATE_SUBLEVELS and INTEGRATE_SUBLEVELS routines test
  ! ======================================================
  WRITE( *,'(/5x,"CREATE_SUBLEVELS and INTEGRATE_SUBLEVELS routines test....")' )
  
  ! Interpolate the data on SubLevels
  ! ---------------------------------
  Error_Status = Create_Sublevels( Pressure,    &  ! Input
                                   Temperature, &  ! Input
                                   Absorber,    &  ! Input
                                   N_PER_LAYER, &  ! Input
                                   Sublevel_p,  &  ! Output
                                   Sublevel_t,  &  ! Output
                                   Sublevel_aa  )  ! Output
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error in Create_Sublevels', &
                          FAILURE )
    STOP
  END IF

  ! Calculate the sublevel geopotential heights
  ! -------------------------------------------
  ! Convert water vapor ppmv to partial pressure
  Sublevel_H2O_pp = PPMV_to_PP( Sublevel_p, Sublevel_aa(:,H2O_J_INDEX) )
  IF ( ANY(Sublevel_H2O_pp < ZERO) ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error converting water vapor ppmv to partial pressure', &
                          FAILURE )
    STOP
  END IF
  ! Calculate heights
  Error_Status = Geopotential_Height( Sublevel_p,          &  ! Input
                                      Sublevel_t,          &  ! Input
                                      Sublevel_H2O_pp,     &  ! Input
                                      Sublevel_z,          &  ! Output
                                      Surface_Height=ZERO, &  ! Optional input
                                      Gravity_Correction=1 )  ! Optional input
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error computing geopotential heights', &
                          FAILURE )
    STOP
  END IF

  ! Integrate the profile data
  ! --------------------------
  Error_Status = Integrate_Sublevels( Sublevel_z,  &  ! Input
                                      Sublevel_p,  &  ! Input
                                      Sublevel_t,  &  ! Input
                                      Sublevel_aa, &  ! Input
                                      N_PER_LAYER, &  ! Input
                                      H2O_J_INDEX, &  ! Input
                                      Layer_p,     &  ! Output
                                      Layer_t,     &  ! Output
                                      Layer_aa     )  ! Output
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error integrating the sublevel profile data', &
                          FAILURE )
    STOP
  END IF

  ! Output the data for comparison
  ! ------------------------------
  ! Output some header info
  WRITE( FileID,'("! ======================================================")' )
  WRITE( FileID,'("! CREATE_SUBLEVELS and INTEGRATE_SUBLEVELS routines test")' )
  WRITE( FileID,'("! ======================================================")' )
  ! Output dimension information
  WRITE( FileID,'("! Dimensions; N_SUBLEVELS, N_LEVELS, N_ABSORBERS",/,3(2x,i5))' ) &
                N_SUBLEVELS, N_LEVELS, N_ABSORBERS
  ! Output a header for the sublevel data
  WRITE( FileID,FMT='( "     Z(km) ",2x,&
                      &"     PF(mb)",3x,&
                      &"     T(K)  ",2x,&
                      &"  H2O(ppmv)",5x,&
                      &"  O3(ppmv) ",3x,&
                      &"  N2O(ppmv)",5x,&
                      &"  CO(ppmv) ",3x,&
                      &"  CH4(ppmv)",5x,&
                      &"  N2(ppmv)", &
                      &/, 128( "-" ) )', &
                IOSTAT=IO_Status )
  IF ( IO_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error writing interpolation results header to output file '//FILENAME, &
                          FAILURE )
    STOP
  END IF
  ! Output the sublevel data
  DO k = 1, N_SUBLEVELS
    WRITE( FileID,FMT='(2x,f8.3,4x,f12.7,4x,f7.3,6(2x,es13.6))', &
                  IOSTAT=IO_Status ) Sublevel_z(k)/1000.0_fp, &
                                     Sublevel_p(k), &
                                     Sublevel_t(k), &
                                     ( Sublevel_aa(k,j), j=1,N_ABSORBERS )
    IF ( IO_Status /= 0 ) THEN
      WRITE( Message,'("Error writing interpolation results for sublevel #",i0,&
                      &" to output file ",a)' ) k, FILENAME
      CALL Display_Message( PROGRAM_NAME, &
                            TRIM(Message), &
                            FAILURE )
      STOP
    END IF
  END DO
  ! Output a header for the integration results
  WRITE( FileID,FMT='(2/,"     PL(mb)  ",&
                        &"     TL(K)     ",&
                        &"  H2OL(kmol/cm2) ",&
                        &"  O3L(kmol/cm2)",&
                        &"  N2OL(kmol/cm2) ",&
                        &"  COL(kmol/cm2)",&
                        &"  CH4L(kmol/cm2) ",&
                        &"  N2L(kmol/cm2)", &
                        &/,124( "-" ) )', &
                IOSTAT=IO_Status )
  IF ( IO_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error writing integration results header to output file '//FILENAME, &
                          FAILURE )
    STOP
  END IF
  ! Output the layer integrated data
  DO k = 1, N_LAYERS
    WRITE( FileID,FMT   ='(2x,f10.5,3x,2x,f7.3,3x,6(3x,es13.6))', &
                  IOSTAT=IO_Status ) Layer_p(k), &
                                     Layer_t(k), &
                                     ( Layer_aa(k,j), j=1,N_ABSORBERS )
    IF ( IO_Status /= 0 ) THEN
      WRITE( Message,'("Error writing integration results for layer #",i0,&
                      &" to output file ",a)' ) k, FILENAME
      CALL Display_Message( PROGRAM_NAME, &
                            TRIM(Message), &
                            FAILURE )
      STOP
    END IF
  END DO


  ! ===============================
  ! EFFECTIVE_LAYER_TP routine test
  ! ===============================
  WRITE( *,'(/5x,"EFFECTIVE_LAYER_TP routine test....")' )
  
  ! Convert the water vapour units to partial pressure
  ! --------------------------------------------------
  H2O_pp = PPMV_to_PP( Pressure, Absorber(:,H2O_J_INDEX) )
  IF ( ANY(H2O_pp < ZERO) ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error converting water vapor ppmv to partial pressure '//&
                          'for EFFECTIVE_LAYER_TP routine test', &
                          FAILURE )
    STOP
  END IF
  
  ! Compute the geopotential heights
  ! --------------------------------
  Error_Status = Geopotential_Height( Pressure, Temperature, H2O_pp, z )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error calculating the geopotential heights '//&
                          'for EFFECTIVE_LAYER_TP routine test', &
                          FAILURE )
    STOP
  END IF
  
  ! Compute the effective layer pressure and temperature
  ! ----------------------------------------------------
  Error_Status = Effective_Layer_TP( z, Pressure, Temperature, H2O_pp, &
                                     Effective_Layer_p, Effective_Layer_t )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error calculating the effective layer pressure and temperature', &
                          FAILURE )
    STOP
  END IF
  
  ! Output results
  ! --------------
  ! Output some header info
  WRITE( FileID,'(3/)' )
  WRITE( FileID,'("! ===============================")' )
  WRITE( FileID,'("! EFFECTIVE_LAYER_TP routine test")' )
  WRITE( FileID,'("! ===============================")' )
  WRITE( FileID,'(/6x,"P(Bot)",6x,"P(Top)",6x,"P(Eff)",8x,"PL",&
                     &8x,"T(Bot)",6x,"T(Top)",6x,"T(Eff)",8x,"TL")' )
  ! Output the effective layer test results
  DO k = 1, N_LAYERS
    WRITE( FileID,'( 2x,8( f12.6 ))' ) &
                  Pressure(k), Pressure(k+1), Effective_Layer_p(k), Layer_p(k), &
                  Temperature(k), Temperature(k+1), Effective_Layer_t(k), Layer_t(k)
  END DO
  
  ! Close the output file
  ! ---------------------
  CLOSE( FileID )

END PROGRAM Level_Layer_Conversion_Test
