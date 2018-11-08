!
! Atmospheric_Properties_Test
!
! Program to test the routines in the Atmospheric_Properties module
! of the Profile_Utility library.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 22-Nov-2004
!                       paul.vandelst@ssec.wisc.edu
!

PROGRAM Atmospheric_Properties_Test

  ! ----------------
  ! Enviroment setup
  ! ----------------
  ! Module usage
  USE Profile_Utility
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'Atmospheric_Properties_Test'
  CHARACTER(*), PARAMETER :: PROGRAM_RCS_ID = &
    '$Id: Atmospheric_Properties_Test.f90 1960 2008-03-18 17:26:10Z paul.vandelst@noaa.gov $'
  ! The dimension
  INTEGER, PARAMETER :: N_LEVELS = 50
  ! The minimum pressure to test the saturation mixing ratio routine
  REAL(fp), PARAMETER :: MIN_PRESSURE = 50.0_fp
  ! The profile data. Absorber amounts are in ppmv.
  ! Original UMBC profile #1 (Tropical Atm)
  REAL(fp), PARAMETER :: PRESSURE(N_LEVELS) = &
  (/ 1.01300e+03_fp,  9.04000e+02_fp,  8.05000e+02_fp,  7.15000e+02_fp, &
     6.33000e+02_fp,  5.59000e+02_fp,  4.92000e+02_fp,  4.32000e+02_fp, &
     3.78000e+02_fp,  3.29000e+02_fp,  2.86000e+02_fp,  2.47000e+02_fp, &
     2.13000e+02_fp,  1.82000e+02_fp,  1.56000e+02_fp,  1.32000e+02_fp, &
     1.11000e+02_fp,  9.37000e+01_fp,  7.89000e+01_fp,  6.66000e+01_fp, &
     5.65000e+01_fp,  4.80000e+01_fp,  4.09000e+01_fp,  3.50000e+01_fp, &
     3.00000e+01_fp,  2.57000e+01_fp,  1.76300e+01_fp,  1.22000e+01_fp, &
     8.52000e+00_fp,  6.00000e+00_fp,  4.26000e+00_fp,  3.05000e+00_fp, &
     2.20000e+00_fp,  1.59000e+00_fp,  1.16000e+00_fp,  8.54000e-01_fp, &
     4.56000e-01_fp,  2.39000e-01_fp,  1.21000e-01_fp,  5.80000e-02_fp, &
     2.60000e-02_fp,  1.10000e-02_fp,  4.40000e-03_fp,  1.72000e-03_fp, &
     6.88000e-04_fp,  2.89000e-04_fp,  1.30000e-04_fp,  6.47000e-05_fp, &
     3.60000e-05_fp,  2.25000e-05_fp /)
  REAL(fp), PARAMETER :: TEMPERATURE(N_LEVELS) = &
  (/ 299.70_fp, 293.70_fp, 287.70_fp, 283.70_fp, &
     277.00_fp, 270.30_fp, 263.60_fp, 257.00_fp, &
     250.30_fp, 243.60_fp, 237.00_fp, 230.10_fp, &
     223.60_fp, 217.00_fp, 210.30_fp, 203.70_fp, &
     197.00_fp, 194.80_fp, 198.80_fp, 202.70_fp, &
     206.70_fp, 210.70_fp, 214.60_fp, 217.00_fp, &
     219.20_fp, 221.40_fp, 227.00_fp, 232.30_fp, &
     237.70_fp, 243.10_fp, 248.50_fp, 254.00_fp, &
     259.40_fp, 264.80_fp, 269.60_fp, 270.20_fp, &
     263.40_fp, 253.10_fp, 236.00_fp, 218.90_fp, &
     201.80_fp, 184.80_fp, 177.10_fp, 177.00_fp, &
     184.30_fp, 190.70_fp, 212.00_fp, 241.60_fp, &
     299.70_fp, 380.00_fp /)
  REAL(fp), PARAMETER :: WATER_VAPOR(N_LEVELS) = &
  (/ 2.593e+04_fp,  1.949e+04_fp,  1.534e+04_fp,  8.600e+03_fp, &
     4.441e+03_fp,  3.346e+03_fp,  2.101e+03_fp,  1.289e+03_fp, &
     7.637e+02_fp,  4.098e+02_fp,  1.912e+02_fp,  7.306e+01_fp, &
     2.905e+01_fp,  9.900e+00_fp,  6.220e+00_fp,  4.000e+00_fp, &
     3.000e+00_fp,  2.900e+00_fp,  2.750e+00_fp,  2.600e+00_fp, &
     2.600e+00_fp,  2.650e+00_fp,  2.800e+00_fp,  2.900e+00_fp, &
     3.200e+00_fp,  3.250e+00_fp,  3.600e+00_fp,  4.000e+00_fp, &
     4.300e+00_fp,  4.600e+00_fp,  4.900e+00_fp,  5.200e+00_fp, &
     5.500e+00_fp,  5.700e+00_fp,  5.900e+00_fp,  6.000e+00_fp, &
     6.000e+00_fp,  6.000e+00_fp,  5.400e+00_fp,  4.500e+00_fp, &
     3.300e+00_fp,  2.100e+00_fp,  1.300e+00_fp,  8.500e-01_fp, &
     5.400e-01_fp,  4.000e-01_fp,  3.400e-01_fp,  2.800e-01_fp, &
     2.400e-01_fp,  2.000e-01_fp /)
  REAL(fp), PARAMETER :: CARBON_DIOXIDE(N_LEVELS) = &
  (/ 3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp, &
     3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp, &
     3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp, &
     3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp, &
     3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp, &
     3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp, &
     3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp, &
     3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp, &
     3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp, &
     3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp,  3.700e+02_fp, &
     3.700e+02_fp,  3.678e+02_fp,  3.588e+02_fp,  3.476e+02_fp, &
     3.027e+02_fp,  2.186e+02_fp,  1.233e+02_fp,  6.727e+01_fp, &
     4.485e+01_fp,  3.924e+01_fp /)
  REAL(fp), PARAMETER :: OZONE(N_LEVELS) = &
  (/ 2.869e-02_fp,  3.150e-02_fp,  3.342e-02_fp,  3.504e-02_fp, &
     3.561e-02_fp,  3.767e-02_fp,  3.989e-02_fp,  4.223e-02_fp, &
     4.471e-02_fp,  5.000e-02_fp,  5.595e-02_fp,  6.613e-02_fp, &
     7.815e-02_fp,  9.289e-02_fp,  1.050e-01_fp,  1.256e-01_fp, &
     1.444e-01_fp,  2.500e-01_fp,  5.000e-01_fp,  9.500e-01_fp, &
     1.400e+00_fp,  1.800e+00_fp,  2.400e+00_fp,  3.400e+00_fp, &
     4.300e+00_fp,  5.400e+00_fp,  7.800e+00_fp,  9.300e+00_fp, &
     9.850e+00_fp,  9.700e+00_fp,  8.800e+00_fp,  7.500e+00_fp, &
     5.900e+00_fp,  4.500e+00_fp,  3.450e+00_fp,  2.800e+00_fp, &
     1.800e+00_fp,  1.100e+00_fp,  6.500e-01_fp,  3.000e-01_fp, &
     1.800e-01_fp,  3.300e-01_fp,  5.000e-01_fp,  5.200e-01_fp, &
     5.000e-01_fp,  4.000e-01_fp,  2.000e-01_fp,  5.000e-02_fp, &
     5.000e-03_fp,  5.000e-04_fp /)
  REAL(fp), PARAMETER :: CARBON_MONOXIDE(N_LEVELS) = &
  (/ 1.441e-01_fp,  1.301e-01_fp,  1.205e-01_fp,  1.185e-01_fp, &
     1.300e-01_fp,  1.433e-01_fp,  1.394e-01_fp,  1.467e-01_fp, &
     9.681e-02_fp,  4.414e-02_fp,  4.255e-02_fp,  4.168e-02_fp, &
     3.936e-02_fp,  4.099e-02_fp,  3.844e-02_fp,  3.595e-02_fp, &
     3.398e-02_fp,  3.555e-02_fp,  2.084e-02_fp,  2.102e-02_fp, &
     2.016e-02_fp,  1.808e-02_fp,  1.507e-02_fp,  1.872e-02_fp, &
     1.781e-02_fp,  2.180e-02_fp,  1.975e-02_fp,  2.079e-02_fp, &
     1.893e-02_fp,  2.145e-02_fp,  2.706e-02_fp,  3.170e-02_fp, &
     3.084e-02_fp,  4.384e-02_fp,  5.032e-02_fp,  6.413e-02_fp, &
     9.016e-02_fp,  1.318e-01_fp,  2.354e-01_fp,  3.243e-01_fp, &
     7.947e-01_fp,  1.873e+00_fp,  3.739e+00_fp,  7.558e+00_fp, &
     1.412e+01_fp,  1.841e+01_fp,  2.602e+01_fp,  3.510e+01_fp, &
     4.369e+01_fp,  5.532e+01_fp /)
  REAL(fp), PARAMETER :: METHANE(N_LEVELS) = &
  (/ 1.800e+00_fp,  1.800e+00_fp,  1.800e+00_fp,  1.800e+00_fp, &
     1.800e+00_fp,  1.800e+00_fp,  1.800e+00_fp,  1.799e+00_fp, &
     1.797e+00_fp,  1.793e+00_fp,  1.784e+00_fp,  1.774e+00_fp, &
     1.760e+00_fp,  1.742e+00_fp,  1.722e+00_fp,  1.699e+00_fp, &
     1.675e+00_fp,  1.644e+00_fp,  1.610e+00_fp,  1.567e+00_fp, &
     1.508e+00_fp,  1.435e+00_fp,  1.347e+00_fp,  1.261e+00_fp, &
     1.184e+00_fp,  1.117e+00_fp,  1.045e+00_fp,  9.673e-01_fp, &
     8.788e-01_fp,  7.899e-01_fp,  7.007e-01_fp,  5.970e-01_fp, &
     4.885e-01_fp,  3.845e-01_fp,  2.936e-01_fp,  2.224e-01_fp, &
     1.748e-01_fp,  1.588e-01_fp,  1.588e-01_fp,  1.588e-01_fp, &
     1.588e-01_fp,  1.588e-01_fp,  1.588e-01_fp,  1.482e-01_fp, &
     1.376e-01_fp,  1.271e-01_fp,  1.165e-01_fp,  1.006e-01_fp, &
     6.353e-02_fp,  3.176e-02_fp /)


  ! ---------
  ! Variables
  ! ---------
  REAL(fp) :: P_Test 
  REAL(fp) :: ppH2O_Test 
  REAL(fp) :: T_Test, T_Answer
  REAL(fp) :: MW_Test, MW_Answer
  REAL(fp) :: Rho_Test
  REAL(fp) :: W_Test
  REAL(fp) :: Tv_Test, Tv_Answer
  REAL(fp) :: Tvi_Test, Tvi_Answer
  REAL(fp) :: Theta_Test, Theta_Answer
  REAL(fp) :: Numerator
  REAL(fp) :: Denominator

  ! Output program header
  ! ---------------------
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to test the routines in the Atmospheric_Properties module '//&
                        'of the Profile_Utility library.', &
                        '$Revision: 1960 $' )

  ! Compute molecular weights
  ! -------------------------
  ppH2O_Test = ONE / ( MW_DRYAIR - MW_H2O )
  P_Test     = 1000.0_fp
  MW_Test    = MW_Air( P_Test, ppH2O_Test )
  MW_Answer  = MW_DRYAIR - 0.001_fp
  CALL Test_Results( MW_Test, MW_Answer, 'MW_Air                       ' )

  ! Compute densities
  ! -----------------
  P_Test  = 1000.0_fp
  T_Test  = 100.0_fp
  MW_Test = R0
  Rho_Test = Density( P_Test, T_Test, MW_Test )
  CALL Test_Results( Rho_Test, ONE, 'Density                      ' )

  ! Compute virtual temperatures
  ! ----------------------------
  T_Test = T0
  W_Test = ONE
  Tv_Test = Virtual_Temperature( T_Test, W_Test )
  Numerator   = EPS + ( G_TO_KG * W_Test )
  Denominator = EPS * ( ONE + ( G_TO_KG * W_Test ) )
  Tv_Answer = T_Test * Numerator / Denominator
  CALL Test_Results( Tv_Test, Tv_Answer, 'Virtual_Temperature   FORWARD' )
  Tvi_Test   = Virtual_Temperature( Tv_Test, W_Test, Inverse = 1 )
  Tvi_Answer = T_Test
  CALL Test_Results( Tvi_Test, Tvi_Answer, 'Virtual_Temperature   INVERSE' )

  ! Compute potential temperatures
  ! ------------------------------
  P_Test = P0
  T_Test = T0
  Theta_Test = Potential_Temperature( T_Test, P_Test )
  Theta_Answer = T0
  CALL Test_Results( Theta_Test, Theta_Answer, 'Potential_Temperature FORWARD' )
  T_Test   = Potential_Temperature( Theta_Test, P_Test, INVERSE = 1 )
  T_Answer = T0
  CALL Test_Results( T_Test, T_Answer, 'Potential_Temperature INVERSE' )


CONTAINS


  ! Simple subroutine to compare data
  ! ---------------------------------
  SUBROUTINE Test_Results( x, y, Comparison )
    REAL(fp),     INTENT(IN) :: x, y
    CHARACTER(*), INTENT(IN) :: Comparison

    INTEGER, PARAMETER ::     MAX_ULPS = 10
    INTEGER, PARAMETER :: WARNING_ULPS = 5

    CHARACTER(256) :: Message
    REAL(fp) :: Value
    REAL(fp) :: Diff
    REAL(fp) :: Diffpc
    INTEGER  :: n, ULP
    LOGICAL  :: Equal
    INTEGER  :: Error_Status

    Diff   = ABS(x-y)
    Value  = MAX(ABS(x),ABS(y))
    Diffpc = 100.0_fp * Diff / Value

    Equal = .FALSE.
    ULP_Loop: DO n = 1, MAX_ULPS
      IF ( Compare_Float( x, y, ULP = n ) ) THEN
        Equal = .TRUE.
        ULP   = n
        EXIT ULP_Loop
      END IF
    END DO ULP_Loop

    IF ( Equal ) THEN
      IF ( ULP <= WARNING_ULPS ) THEN
        Error_Status = INFORMATION
      ELSE
        Error_Status = WARNING
      END IF

      WRITE( Message,'(a," consistent. ULP = ",i2,", %|Diff| = ",es13.6)' ) &
                     Comparison, ULP, Diffpc
      CALL Display_Message( PROGRAM_NAME, &
                            TRIM(Message), &
                            Error_Status )
    ELSE
      WRITE( Message,'( a," inconsistent.  %|Diff| = ",es13.6)' ) &
                     Comparison, Diffpc
      CALL Display_Message( PROGRAM_NAME, &
                            TRIM(Message), &
                            FAILURE )
      CALL Output_Difference( x, y, Comparison )
    END IF
  END SUBROUTINE Test_Results


  ! Simple subroutine to output data when comparison fails
  SUBROUTINE Output_Difference( New, Old, Units )
    REAL(fp),     INTENT(IN) :: New, Old
    CHARACTER(*), INTENT(IN) :: Units

    REAL(fp) :: Difference
    REAL(fp) :: Percentage_Difference

    Difference = New - Old
    IF ( New > ZERO ) THEN
      Percentage_Difference = 100.0_fp * Difference / New
    ELSE
      Percentage_Difference = 999.0_fp
    END IF

    WRITE( *,'(/5x, "DATA UNITS COMPARISON : ",a)' ) Units
    WRITE( *,'("       NEW          ORIGINAL        Diff.         % Diff")' )
    WRITE( *,'(2x,56("-"))' )
    WRITE( *,'(3(2x,es13.6),2x,f10.5)' ) &
             New, Old, &
             Difference, Percentage_Difference
  END SUBROUTINE Output_Difference

END PROGRAM Atmospheric_Properties_Test
