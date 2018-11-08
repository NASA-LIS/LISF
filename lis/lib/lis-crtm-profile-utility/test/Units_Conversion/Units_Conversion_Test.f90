!
! Units_Conversion_Test
!
! Program to test the routines in the Units_Conversion module of the 
! Profle_Utility library.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 22-Nov-2004
!                       paul.vandelst@ssec.wisc.edu
!

PROGRAM Units_Conversion_Test

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
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'Units_Conversion_Test'
  CHARACTER(*), PARAMETER :: PROGRAM_RCS_ID = &
    '$Id: Units_Conversion_Test.f90 1954 2008-03-17 20:40:21Z paul.vandelst@noaa.gov $'
  ! The dimensions
  INTEGER, PARAMETER :: N_LEVELS = 50
  INTEGER, PARAMETER :: N_LAYERS = N_LEVELS-1
  ! The number of absorbers and their IDs
  INTEGER, PARAMETER :: N_ABSORBERS = 5
  INTEGER, PARAMETER :: MOLECULE_ID(N_ABSORBERS) = (/ ID_H2O, &
                                                      ID_CO2, &
                                                      ID_O3,  &
                                                      ID_CO,  &
                                                      ID_CH4 /)
  CHARACTER(*), PARAMETER :: MOLECULE_NAME(N_ABSORBERS) = (/ 'WATER_VAPOR    ', &
                                                             'CARBON_DIOXIDE ', &
                                                             'OZONE          ', &
                                                             'CARBON_MONOXIDE', &
                                                             'METHANE        ' /)
  ! The minimum pressure to test the relative humidity routines
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
  INTEGER :: Error_Status
  INTEGER :: j, k

  REAL(fp), DIMENSION(N_LEVELS) :: H2O_Pressure
  REAL(fp), DIMENSION(N_LEVELS) :: H2O_ppmv
  REAL(fp), DIMENSION(N_LEVELS) :: H2O_Number_Density
  REAL(fp), DIMENSION(N_LEVELS) :: H2O_Mixing_Ratio
  REAL(fp), DIMENSION(N_LEVELS) :: H2O_Specific_Amount
  REAL(fp), DIMENSION(N_LAYERS) :: H2O_Layer_ppmv
  REAL(fp), DIMENSION(N_LAYERS) :: H2O_kmolcm2
  REAL(fp), DIMENSION(N_LEVELS) :: Height
  REAL(fp), DIMENSION(N_LEVELS) :: dPressure
  INTEGER                       :: pIdx
  REAL(fp), DIMENSION(N_LEVELS) :: Profile_Data
  REAL(fp), DIMENSION(N_LEVELS) :: Mixing_Ratio, New_Mixing_Ratio
  REAL(fp), DIMENSION(N_LEVELS) :: Specific_Amount
  REAL(fp), DIMENSION(N_LEVELS) :: Relative_Humidity
  REAL(fp), DIMENSION(N_LEVELS) :: Partial_Pressure, New_Partial_Pressure
  REAL(fp), DIMENSION(N_LEVELS) :: ppmv
  REAL(fp), DIMENSION(N_LEVELS) :: Number_Density
  REAL(fp), DIMENSION(N_LEVELS) :: Mass_Density
  REAL(fp), DIMENSION(N_LAYERS) :: Layer_Pressure
  REAL(fp), DIMENSION(N_LAYERS) :: Layer_Temperature
  REAL(fp), DIMENSION(N_LAYERS) :: Layer_Thickness
  REAL(fp), DIMENSION(N_LAYERS) :: Layer_ppmv, New_Layer_ppmv
  REAL(fp), DIMENSION(N_LAYERS) :: Layer_kmolcm2

  ! Output program header
  ! ---------------------
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to test the Profile_Utility Units_Conversion '//&
                        ' module routines.', &
                        '$Revision: 1954 $' )

  ! Gather the various water vapour bits
  ! ------------------------------------
  H2O_ppmv            = WATER_VAPOR
  H2O_Pressure        = PPMV_to_PP( PRESSURE, WATER_VAPOR )
  H2O_Number_Density  = PPMV_to_ND( PRESSURE, TEMPERATURE, WATER_VAPOR )
  H2O_Mixing_Ratio    = PPMV_to_MR( WATER_VAPOR )
  H2O_Specific_Amount = MR_to_SA( H2O_Mixing_Ratio )

  ! Get set up for the kmol/cm^2 conversions
  ! ----------------------------------------
  ! Compute the geopotential heights
  Error_Status = Geopotential_Height( PRESSURE, &
                                      TEMPERATURE, &
                                      H2O_Pressure, &
                                      Height )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error calculating the geopotential heights', &
                          FAILURE )
    STOP
  END IF
  ! Compute the layer thickness
  Layer_Thickness = Height( 2:N_LEVELS ) - Height( 1:N_LEVELS-1 )
  ! Compute the effective layer pressure and temperature
  Error_Status = Effective_Layer_TP( Height, &
                                     PRESSURE, &
                                     TEMPERATURE, &
                                     H2O_Pressure, &
                                     Layer_Pressure, &
                                     Layer_Temperature )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error calculating the effective layer pressure and temperature', &
                          FAILURE )
    STOP
  END IF
  ! Compute the water vapor kmol/cm2
  H2O_Layer_ppmv = 0.5_fp * ( H2O_ppmv(1:N_LEVELS-1) + H2O_ppmv(2:N_LEVELS) )
  H2O_kmolcm2 = PPMV_to_KMOL( Layer_Pressure, &
                              Layer_Temperature, &
                              Layer_Thickness, &
                              H2O_Layer_ppmv )


  ! Begin loop over absorber
  ! ------------------------
  Absorber_Loop: DO j = 1, N_ABSORBERS

    ! Load an absorber profile
    ! ------------------------
    SELECT CASE ( MOLECULE_ID(j) )
      CASE ( ID_H2O )
        Profile_Data = WATER_VAPOR
      CASE ( ID_CO2 )
        Profile_Data = CARBON_DIOXIDE
      CASE ( ID_O3 )
        Profile_Data = OZONE
      CASE ( ID_CO )
        Profile_Data = CARBON_MONOXIDE
      CASE ( ID_CH4 )
        Profile_Data = METHANE
    END SELECT
    WRITE( *,'(//5x,"Molecule: ",a)' ) MOLECULE_NAME(j)


    ! Test routine pairs for consistency
    ! ----------------------------------
    WRITE( *,'(/10x,"Testing conversion routines for consistency...")' )

    ! ppmv <-> partial pressure
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      Partial_Pressure = PPMV_to_PP( PRESSURE, Profile_Data )
    ELSE
      Partial_Pressure = PPMV_to_PP( PRESSURE, Profile_Data, Water_Vapor=H2O_Pressure )
    END IF
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      ppmv = PP_to_PPMV( PRESSURE, Partial_Pressure )
    ELSE
      ppmv = PP_to_PPMV( PRESSURE, Partial_Pressure, Water_Vapor=H2O_ppmv )
    END IF
    CALL Test_Results( ppmv, Profile_Data, 'PP<->PPMV  ' )

    ! ppmv <-> mixing ratio
    Mixing_Ratio = PPMV_to_MR( Profile_Data, Molecule_ID=MOLECULE_ID(j) )
    ppmv = MR_to_PPMV( Mixing_Ratio, Molecule_ID=MOLECULE_ID(j) )
    CALL Test_Results( ppmv, Profile_Data, 'MR<->PPMV  ' )

    ! mixing ratio <-> specific amounts
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      Specific_Amount = MR_to_SA( Mixing_Ratio )
    ELSE
      Specific_Amount = MR_to_SA( Mixing_Ratio, Water_Vapor=H2O_Specific_Amount )
    END IF
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      New_Mixing_Ratio = SA_to_MR( Specific_Amount )
    ELSE
      New_Mixing_Ratio = SA_to_MR( Specific_Amount, Water_Vapor=H2O_Mixing_Ratio )
    END IF
    CALL Test_Results( New_Mixing_Ratio, Mixing_Ratio, 'MR<->SA    ' )

    ! Partial pressure <-> mixing ratio
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      Partial_Pressure = MR_to_PP( PRESSURE, Mixing_Ratio )
    ELSE
      Partial_Pressure = MR_to_PP( PRESSURE, Mixing_Ratio, &
                                   Molecule_ID=MOLECULE_ID(j), &
                                   Water_Vapor=H2O_Pressure )
    END IF
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      New_Mixing_Ratio = PP_to_MR( PRESSURE, Partial_Pressure )
    ELSE
      New_Mixing_Ratio = PP_to_MR( PRESSURE, Partial_Pressure, &
                                   Molecule_ID=MOLECULE_ID(j), &
                                   Water_Vapor=H2O_Mixing_Ratio )
    END IF
    CALL Test_Results( New_Mixing_Ratio, Mixing_Ratio, 'MR<->PP    ' )

    ! Mixing ratio <-> relative humidity (Water vapor only)
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      Relative_Humidity = MR_to_RH( PRESSURE, TEMPERATURE, Mixing_Ratio, &
                                    Min_Pressure=MIN_PRESSURE )
      New_Mixing_Ratio  = RH_to_MR( PRESSURE, TEMPERATURE, Relative_Humidity, &
                                    Min_Pressure=MIN_PRESSURE )
      dPressure = PRESSURE - MIN_PRESSURE
      pIdx = MINLOC(dPressure, MASK=dPressure > ZERO, DIM=1)
      CALL Test_Results( New_Mixing_Ratio(1:pIdx), Mixing_Ratio(1:pIdx), 'MR<->RH    ' )
    END IF

    ! Partial pressure <-> mass density
    Mass_Density = PP_to_MD( Partial_Pressure, TEMPERATURE, &
                             Molecule_ID=MOLECULE_ID(j) )
    New_Partial_Pressure = MD_to_PP( Mass_Density, TEMPERATURE, &
                                    Molecule_ID=MOLECULE_ID(j) ) 
    CALL Test_Results( New_Partial_Pressure, Partial_Pressure, 'PP<->MD    ' )

    ! Partial pressure <-> number density
    Number_Density = PP_to_ND( Partial_Pressure, TEMPERATURE )
    New_Partial_Pressure = ND_to_PP( Number_Density, TEMPERATURE ) 
    CALL Test_Results( New_Partial_Pressure, Partial_Pressure, 'PP<->ND    ' )


    ! ppmv <-> number density
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      Number_Density = PPMV_to_ND( PRESSURE, TEMPERATURE, Profile_Data )
      ppmv = ND_to_PPMV( PRESSURE, TEMPERATURE, Number_Density )
    ELSE
      Number_Density = PPMV_to_ND( PRESSURE, Temperature, Profile_Data, &
                                   Water_Vapor=H2O_Number_Density )
      ppmv = ND_to_PPMV( PRESSURE, TEMPERATURE, Number_Density, &
                         Water_Vapor=H2O_ppmv )
    END IF
    CALL Test_Results( ppmv, Profile_Data, 'PPMV<->ND  ' )

    ! ppmv <-> kmol_per_cm2
    Layer_ppmv = 0.5_fp * ( Profile_Data(1:N_LEVELS-1) + Profile_Data(2:N_LEVELS) )
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      Layer_kmolcm2 = PPMV_to_KMOL( Layer_Pressure, &
                                    Layer_Temperature, &
                                    Layer_Thickness, &
                                    Layer_ppmv )
      New_Layer_ppmv = KMOL_to_PPMV( Layer_Pressure, &
                                     Layer_Temperature, &
                                     Layer_Thickness, &
                                     Layer_kmolcm2 )
    ELSE
      Layer_kmolcm2 = PPMV_to_KMOL( Layer_Pressure, &
                                    Layer_Temperature, &
                                    Layer_Thickness, &
                                    Layer_ppmv, &
                                    Water_Vapor = H2O_kmolcm2 )
      New_Layer_ppmv = KMOL_to_PPMV( Layer_Pressure, &
                                     Layer_Temperature, &
                                     Layer_Thickness, &
                                     Layer_kmolcm2, &
                                     Water_Vapor = H2O_Layer_ppmv )
    END IF
    CALL Test_Results( New_Layer_ppmv, Layer_ppmv, 'PPMV<->KMOL' )


    ! Convert data units in a chain
    ! -----------------------------
    WRITE( *,'(/10x,"Testing conversion routines in a chain....")' )
    WRITE( *,'(15x,"Converting PPMV -> MR....")' )
    Mixing_Ratio = PPMV_to_MR( Profile_Data, Molecule_ID=MOLECULE_ID(j) )
    
    WRITE( *,'(15x,"Converting MR -> PP....")' )
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      Partial_Pressure = MR_to_PP( PRESSURE, Mixing_Ratio )
    ELSE
      Partial_Pressure = MR_to_PP( PRESSURE, Mixing_Ratio, &
                                   Molecule_ID=MOLECULE_ID(j), &
                                   Water_Vapor=H2O_Pressure )
    END IF
    
    WRITE( *,'(15x,"Converting PP -> ND....")' )
    Number_Density = PP_to_ND( Partial_Pressure, TEMPERATURE )

    WRITE( *,'(15x,"Converting ND -> PPMV....")' )
    IF ( MOLECULE_ID(j) == ID_H2O ) THEN
      ppmv = ND_to_PPMV( PRESSURE, TEMPERATURE, &
                         Number_Density )
    ELSE
      ppmv = ND_to_PPMV( PRESSURE, TEMPERATURE, &
                         Number_Density, &
                         Water_Vapor=H2O_ppmv )
    END IF

    CALL Test_Results( ppmv, Profile_Data, 'PPMV->MR->PP->ND->PPMV' )

  END DO Absorber_Loop


CONTAINS


  ! Simple subroutine to compare data
  ! ---------------------------------
  SUBROUTINE Test_Results( x, y, Comparison )
    REAL(fp),     INTENT(IN) :: x(:), y(:)
    CHARACTER(*), INTENT(IN) :: Comparison

    INTEGER, PARAMETER ::     MAX_ULPS = 10
    INTEGER, PARAMETER :: WARNING_ULPS = 5

    CHARACTER(256) :: Message
    INTEGER  :: MaxIdx
    REAL(fp) :: MaxValue
    REAL(fp) :: MaxDiff
    REAL(fp) :: MaxDiffpc
    INTEGER  :: n, ULP
    LOGICAL  :: Equal
    INTEGER  :: Error_Status

    MaxDiff   = MAXVAL(ABS(x-y))
    MaxIdx    = MAXLOC(ABS(x-y), DIM=1)
    MaxValue  = MAX(ABS(x(MaxIdx)), ABS(y(MaxIdx)))
    MaxDiffpc = 100.0_fp * MaxDiff / MaxValue

    Equal = .FALSE.
    ULP_Loop: DO n = 1, MAX_ULPS
      IF ( ALL( Compare_Float( x, y, ULP = n ) ) ) THEN
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

      WRITE( Message,'(a," consistent. ULP = ",i2,", %MAX(|Diff|) = ",es13.6)' ) &
                     Comparison, ULP, MaxDiffpc
      CALL Display_Message( PROGRAM_NAME, &
                            TRIM( Message ), &
                            Error_Status )
    ELSE
      WRITE( Message,'( a," inconsistent.  %MAX(|Diff|) = ",es13.6)' ) &
                     Comparison, MaxDiffpc
      CALL Display_Message( PROGRAM_NAME, &
                            TRIM( Message ), &
                            FAILURE )
      CALL Output_Difference( x, y, Comparison )
    END IF
  END SUBROUTINE Test_Results


  ! Simple subroutine to output data when comparison fails
  SUBROUTINE Output_Difference( New, Old, Units )
    REAL(fp),     INTENT(IN) :: New(:), Old(:)
    CHARACTER(*), INTENT(IN) :: Units

    REAL(fp), DIMENSION(SIZE(New)) :: Difference
    REAL(fp), DIMENSION(SIZE(New)) :: Percentage_Difference

    Difference = New - Old
    WHERE( New > ZERO )
      Percentage_Difference = 100.0_fp * Difference / New
    ELSEWHERE
      Percentage_Difference = 999.0_fp
    END WHERE

    WRITE( *,'(/5x, "DATA UNITS COMPARISON : ",a)' ) Units
    WRITE( *,'("       NEW          ORIGINAL        Diff.         % Diff")' )
    WRITE( *,'(2x,56("-"))' )
    DO k = 1, SIZE(New)
      WRITE( *,'(3(2x,es13.6),2x,f10.5)' ) &
                New(k), Old(k), &
                Difference(k), Percentage_Difference(k)
    END DO
  END SUBROUTINE Output_Difference

END PROGRAM Units_Conversion_Test
