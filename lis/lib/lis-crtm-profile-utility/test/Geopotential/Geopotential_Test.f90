!
! Geopotential_Height_Test
!
! Program to test the routines in the Geopotential module
! of the Profile_Utility library.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 21-Jan-2003
!                       paul.vandelst@ssec.wisc.edu
!

PROGRAM Geopotential_Test

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
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'Geopotential_Test'
  CHARACTER(*), PARAMETER :: PROGRAM_RCS_ID = &
  '$Id: Geopotential_Test.f90 1960 2008-03-18 17:26:10Z paul.vandelst@noaa.gov $'
  ! Output filename
  CHARACTER(*), PARAMETER :: FILENAME = 'Geopotential_Test.asc'
  ! Define dimension
  INTEGER, PARAMETER :: N_LEVELS = 101
  ! Define test profile pressure
  REAL(fp), PARAMETER :: PRESSURE(N_LEVELS) = &
    (/ 1100.000_fp,1070.917_fp,1042.232_fp,1013.948_fp, 986.067_fp, &
        958.591_fp, 931.524_fp, 904.866_fp, 878.620_fp, 852.788_fp, &
        827.371_fp, 802.371_fp, 777.790_fp, 753.628_fp, 729.886_fp, &
        706.565_fp, 683.667_fp, 661.192_fp, 639.140_fp, 617.511_fp, &
        596.306_fp, 575.525_fp, 555.167_fp, 535.232_fp, 515.720_fp, &
        496.630_fp, 477.961_fp, 459.712_fp, 441.882_fp, 424.470_fp, &
        407.474_fp, 390.893_fp, 374.724_fp, 358.966_fp, 343.618_fp, &
        328.675_fp, 314.137_fp, 300.000_fp, 286.262_fp, 272.919_fp, &
        259.969_fp, 247.409_fp, 235.234_fp, 223.441_fp, 212.028_fp, &
        200.989_fp, 190.320_fp, 180.018_fp, 170.078_fp, 160.496_fp, &
        151.266_fp, 142.385_fp, 133.846_fp, 125.646_fp, 117.777_fp, &
        110.237_fp, 103.017_fp,  96.114_fp,  89.520_fp,  83.231_fp, &
         77.240_fp,  71.540_fp,  66.125_fp,  60.990_fp,  56.126_fp, &
         51.528_fp,  47.188_fp,  43.100_fp,  39.257_fp,  35.650_fp, &
         32.274_fp,  29.121_fp,  26.183_fp,  23.453_fp,  20.922_fp, &
         18.585_fp,  16.432_fp,  14.456_fp,  12.649_fp,  11.004_fp, &
          9.512_fp,   8.165_fp,   6.957_fp,   5.878_fp,   4.920_fp, &
          4.077_fp,   3.340_fp,   2.701_fp,   2.153_fp,   1.687_fp, &
          1.297_fp,   0.975_fp,   0.714_fp,   0.506_fp,   0.345_fp, &
          0.224_fp,   0.137_fp,   0.077_fp,   0.038_fp,   0.016_fp, &
          0.005_fp /)
  ! Define test profile temperature
  REAL(fp), PARAMETER :: TEMPERATURE(N_LEVELS) = &
    (/ 304.043_fp,302.630_fp,301.199_fp,299.749_fp,298.280_fp, &
       296.790_fp,295.281_fp,293.750_fp,292.227_fp,290.683_fp, &
       289.118_fp,287.590_fp,286.540_fp,285.475_fp,284.395_fp, &
       283.047_fp,281.235_fp,279.397_fp,277.531_fp,275.665_fp, &
       273.782_fp,271.870_fp,269.939_fp,268.020_fp,266.071_fp, &
       264.092_fp,262.131_fp,260.155_fp,258.148_fp,256.118_fp, &
       254.067_fp,251.983_fp,249.880_fp,247.807_fp,245.698_fp, &
       243.553_fp,241.422_fp,239.252_fp,237.043_fp,234.797_fp, &
       232.509_fp,230.178_fp,227.958_fp,225.700_fp,223.408_fp, &
       221.164_fp,218.876_fp,216.524_fp,214.055_fp,211.535_fp, &
       209.083_fp,206.692_fp,204.249_fp,201.792_fp,199.292_fp, &
       196.910_fp,196.031_fp,195.130_fp,195.862_fp,197.557_fp, &
       199.289_fp,201.053_fp,202.874_fp,204.840_fp,206.863_fp, &
       208.960_fp,211.116_fp,213.323_fp,215.232_fp,216.717_fp, &
       218.157_fp,219.623_fp,221.135_fp,222.759_fp,224.456_fp, &
       226.216_fp,228.013_fp,229.857_fp,231.780_fp,233.852_fp, &
       236.043_fp,238.355_fp,240.821_fp,243.424_fp,246.229_fp, &
       249.223_fp,252.505_fp,256.009_fp,259.759_fp,263.815_fp, &
       267.901_fp,269.940_fp,268.260_fp,264.528_fp,258.953_fp, &
       251.472_fp,239.120_fp,225.489_fp,209.888_fp,192.205_fp, &
       178.174_fp /)
  ! Define test profile water vapor in ppmv
  REAL(fp), PARAMETER :: WATER_VAPOR_PPMV(N_LEVELS) = &
    (/ 30590.990_fp,29075.215_fp,27539.310_fp,25982.915_fp,24405.609_fp, &
       22806.964_fp,21186.668_fp,19544.166_fp,18471.102_fp,17403.377_fp, &
       16320.759_fp,15154.037_fp,13385.207_fp,11591.185_fp, 9771.420_fp, &
        8194.816_fp, 7070.010_fp, 5928.731_fp, 4770.583_fp, 4222.798_fp, &
        3915.026_fp, 3602.601_fp, 3278.904_fp, 2922.299_fp, 2560.159_fp, &
        2192.339_fp, 1920.250_fp, 1677.194_fp, 1430.213_fp, 1219.825_fp, &
        1059.069_fp,  895.642_fp,  741.512_fp,  632.000_fp,  520.614_fp, &
         408.258_fp,  337.651_fp,  265.787_fp,  192.629_fp,  153.473_fp, &
         114.298_fp,   74.393_fp,   58.556_fp,   43.271_fp,   28.493_fp, &
          21.983_fp,   15.342_fp,    9.639_fp,    8.283_fp,    6.898_fp, &
           5.810_fp,    5.006_fp,    4.185_fp,    3.715_fp,    3.342_fp, &
           2.996_fp,    2.956_fp,    2.915_fp,    2.860_fp,    2.797_fp, &
           2.731_fp,    2.663_fp,    2.600_fp,    2.600_fp,    2.602_fp, &
           2.628_fp,    2.666_fp,    2.751_fp,    2.826_fp,    2.888_fp, &
           3.058_fp,    3.210_fp,    3.244_fp,    3.335_fp,    3.441_fp, &
           3.551_fp,    3.676_fp,    3.816_fp,    3.961_fp,    4.086_fp, &
           4.208_fp,    4.336_fp,    4.473_fp,    4.618_fp,    4.774_fp, &
           4.939_fp,    5.118_fp,    5.312_fp,    5.513_fp,    5.664_fp, &
           5.829_fp,    5.957_fp,    6.000_fp,    6.000_fp,    6.000_fp, &
           5.943_fp,    5.509_fp,    4.847_fp,    3.868_fp,    2.623_fp, &
           1.412_fp /)


  ! ---------
  ! Variables
  ! ---------
  INTEGER :: FileID
  INTEGER :: IO_Status
  INTEGER :: Error_Status
  REAL(fp) :: Water_Vapor_Pressure(N_LEVELS)
  REAL(fp) :: Height(N_LEVELS)
  REAL(fp) :: Height_with_Gravity_Correction(N_LEVELS)
  REAL(fp) :: Height_at_Latitude_60(N_LEVELS)

  ! Output program header
  ! ---------------------
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to test the routines in the Geopotential module '//&
                        'of the Profile_Utility library.', &
                        '$Revision: 1960 $' )

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

  ! Convert water vapour units from ppmv to hPa
  ! -------------------------------------------
  Water_Vapor_Pressure = PPMV_to_PP( PRESSURE, WATER_VAPOR_PPMV )
  IF ( ANY(Water_Vapor_Pressure < ZERO) ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error converting water vapor ppmv to partial pressure', &
                          FAILURE )
    STOP
  END IF

  ! Calculate the geopotential heights
  ! ----------------------------------
  ! Calculate the plain old heights
  WRITE( *,'(/5x,"Computing geopotential heights...")' )
  Error_Status = Geopotential_Height( PRESSURE, &
                                      TEMPERATURE, &
                                      Water_Vapor_Pressure, &
                                      Height )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Geopotential height calculation failed', &
                          FAILURE )
    STOP
  END IF
  WRITE( FileID,'(5x,"Z(plain) (m):")' )
  WRITE( FileID,'(5f10.1)' ) Height
  ! Calculate the heights, but with a gravity profile correction
  WRITE( *,'(/5x,"Computing geopotential heights with gravity correction...")' )
  Error_Status = Geopotential_Height( PRESSURE, &
                                      TEMPERATURE, &
                                      Water_Vapor_Pressure, &
                                      Height_with_Gravity_Correction, &
                                      Gravity_Correction = 1 )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Geopotential height calculation with gravity correction failed', &
                          FAILURE )
    STOP
  END IF
  WRITE( FileID,'(/5x,"Z(plain) - Z(gcorr) (m):")' )
  WRITE( FileID,'(5f10.2)' ) Height - Height_with_Gravity_Correction
  ! Calculate the heights, but with the gravity profile
  ! correction and at a latitude of 60deg.
  WRITE( *,'(/5x,"Computing geopotential heights with gravity correction at 60deg. latitude...")' )
  Error_Status = Geopotential_Height( PRESSURE, &
                                      TEMPERATURE, &
                                      Water_Vapor_Pressure, &
                                      Height_at_Latitude_60, &
                                      Gravity_Correction = 1, &
                                      Latitude = 60.0_fp )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Geopotential height calculation @ 60deg. latitude failed', &
                          FAILURE )
    STOP
  END IF
  WRITE( FileID,'(/5x,"Z(plain) - Z(gcorr@60N) (m):")' )
  WRITE( FileID,'(5f10.2)' ) Height - Height_at_Latitude_60

  ! Close the output file
  ! ---------------------
  CLOSE( FileID )

END PROGRAM Geopotential_Test
