!
! Atmospheric_Properties
!
! Module containing utility routines to calculate various and sundry
! atmospheric properties.
!
! CONTAINS:
!       MW_Air:                   Function to calculate the effective, water
!                                 vapor weighted molecular weight of air.
!
!       Density:                  Function to calculate gas density using
!                                 the ideal gas law.
!
!       SVP_Water:                Function to calculate the saturation vapor
!                                 pressure over water.
!
!       SVP_Ice:                  Function to calculate the saturation vapor
!                                 pressure over ice.
!
!       Saturation_Mixing_Ratio:  Function to calculate the saturation mixing
!                                 ratio.
!
!       Virtual_Temperature:      Function to calculate the virtual
!                                 temperature, or its inverse.
!
!       Potential_Temperature:    Function to calculate the potential
!                                 temperature, or its inverse.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 01-May-2000
!                       paul.vandelst@noaa.gov
!

MODULE Atmospheric_Properties

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Modules used
  USE Type_Kinds, ONLY: fp
  USE Message_Handler
  USE Profile_Utility_Parameters
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: MW_Air
  PUBLIC :: Density
  PUBLIC :: SVP_Water
  PUBLIC :: SVP_Ice
  PUBLIC :: Saturation_Mixing_Ratio
  PUBLIC :: Virtual_Temperature
  PUBLIC :: Potential_Temperature


  ! ---------------------
  ! Procedure overloading
  ! ---------------------
  INTERFACE MW_Air
    MODULE PROCEDURE MW_Air_scalar
    MODULE PROCEDURE MW_Air_rank1
  END INTERFACE MW_Air

  INTERFACE Density
    MODULE PROCEDURE Density_scalar
    MODULE PROCEDURE Density_rank1a   ! Scalar molecular weight
    MODULE PROCEDURE Density_rank1b   ! Rank-1 molecular weight
  END INTERFACE Density

  INTERFACE SVP_Water
    MODULE PROCEDURE SVPw_scalar
    MODULE PROCEDURE SVPw_rank1
  END INTERFACE SVP_Water

  INTERFACE SVP_Ice
    MODULE PROCEDURE SVPi_scalar
    MODULE PROCEDURE SVPi_rank1
  END INTERFACE SVP_Ice

  INTERFACE Saturation_Mixing_Ratio
    MODULE PROCEDURE SMR_scalar
    MODULE PROCEDURE SMR_rank1
  END INTERFACE Saturation_Mixing_Ratio

  INTERFACE Virtual_Temperature
    MODULE PROCEDURE Tv_scalar
    MODULE PROCEDURE Tv_rank1
  END INTERFACE Virtual_Temperature

  INTERFACE Potential_Temperature
    MODULE PROCEDURE Theta_scalar
    MODULE PROCEDURE Theta_rank1
  END INTERFACE Potential_Temperature


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_RCS_ID = &
  '$Id: Atmospheric_Properties.f90 2018 2008-04-01 20:12:33Z paul.vandelst@noaa.gov $'
  ! Keyword argument set value
  INTEGER, PARAMETER :: SET = 1


CONTAINS


!--------------------------------------------------------------------------------
!S+
! NAME:
!       MW_Air
!
! PURPOSE:
!       Function to calculate the effective, water vapor weighted molecular
!       weight of air.
!
! CALLING SEQUENCE:
!       Molecular_Weight = MW_Air( Pressure,               &  ! Input
!                                  Water_Vapor_Pressure,   &  ! Input
!                                  Message_Log=Message_Log )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:              Total atmospheric pressure
!                              UNITS:      hectoPascals, hPa
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Scalar or Rank-1
!                              ATTRIBUTES: INTENT(IN)
!
!       Water_Vapor_Pressure:  Water vapor partial pressure
!                              UNITS:      hectoPascals, hPa.
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Same as Pressure
!                              ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:           Character string specifying a filename in which any
!                              Messages will be logged. If not specified, or if an
!                              error occurs opening the log file, the default action
!                              is to output Messages to standard output.
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! FUNCTION RESULT:
!       Molecular_Weight:      The effective molecular weight of air.
!                              If an error occurs, -1.0 is returned.
!                              UNITS:      grams, g
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       The change in the effective molecular weight of dry air
!       due to water vapor is given by:
!
!                      pp(h2o) * ( MW(H2O) - MW(DRY_AIR) )
!         d(MW_Air) = -------------------------------------
!                                   Pressure
!
!       and the final result is given by:
!
!         MW(Air) = MW(DRY_AIR) + d(MW_Air)
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 01-May-2000
!                       paul.vandelst@ssec.wisc.edu
!S-
!--------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION MW_Air_scalar( Pressure,             &  ! Input
                          Water_Vapor_Pressure, &  ! Input
                          Message_Log )         &  ! Error messaging
                        RESULT( MW_Air )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Water_Vapor_Pressure
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: MW_Air
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MW_Air'
    REAL(fp),     PARAMETER :: VOLUME_SCALE_FACTOR = ( CM_TO_M )**3
    ! Local variables
    REAL(fp) :: d_MW_Air
    
    ! Setup
    MW_Air = -ONE
    IF ( Pressure < TOLERANCE .OR. Water_Vapor_Pressure < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input pressures < or = 0.0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate change to air molecular weight due to water vapour
    d_MW_Air = Water_Vapor_Pressure * ( MW_H2O - MW_DRYAIR ) / Pressure

    ! Calculate air molecular weight
    MW_Air = MW_DRYAIR + d_MW_Air

  END FUNCTION MW_Air_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION MW_Air_rank1( Pressure,             &  ! Input
                         Water_Vapor_Pressure, &  ! Input
                         Message_Log )         &  ! Error messaging
                       RESULT( MW_Air )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Water_Vapor_Pressure(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: MW_Air(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MW_Air'
    ! Local variables
    INTEGER :: i, n
    
    ! Setup
    MW_Air = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Water_Vapor_Pressure) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate air molecular weight
    DO i = 1, n
      MW_Air( i ) = MW_Air_scalar( Pressure( i ), &
                                   Water_Vapor_Pressure( i ), &
                                   Message_Log=Message_Log )
      IF ( MW_Air( i ) < ZERO ) RETURN
    END DO

  END FUNCTION MW_Air_rank1


!--------------------------------------------------------------------------------
!S+
! NAME:
!       Density
!
! PURPOSE:
!       Function to calculate gas density using the ideal gas law.
!
! CALLING SEQUENCE:
!       Gas_Density = Density( Pressure,               &  ! Input
!                              Temperature,            &  ! Input
!                              Molecular_Weight,       &  ! Input
!                              Message_Log=Message_Log )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Pressure of gas
!                          UNITS:      hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1
!                          ATTRIBUTES: INTENT(IN)
!
!       Temperature:       Temperature of gas
!                          UNITS:      Kelvin
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!
!       Molecular_Weight:  Molecular weight of the gas.
!                          UNITS:      g.mol^-1
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:       Character string specifying a filename in which any
!                          Messages will be logged. If not specified, or if an
!                          error occurs opening the log file, the default action
!                          is to output Messages to standard output.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! FUNCTION RESULT:
!       Gas_Density:       The gas density for the specified conditions.
!                          If an error occurs, -1.0 is returned.
!                          UNITS:      kg.m^-3
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       The density is calculated using the ideal gas equation
!
!                    p . MW
!         Density = --------
!                    R0 . T
!
!       where R0 = universal gas constant.
!
!       Units:
!       ------
!       Pressure           : hPa == 100 Pa == 100 kg.m^-1.s^-2
!       Molecular_Weight   : g.mol^-1 == 0.001 kg.mol^-1
!       MOLAR_GAS_CONSTANT : J.K^-1.mol^-1 == kg.m^2.s^-2.K^-1.mol^-1
!       Temperature        : K
!
!                  100 kg.m^-1.s^-2 . 0.001 kg.mol^-1
!       Density = -----------------------------------
!                     kg.m^2.s^-2.K^-1.mol^-1 . K
!
!                  0.1 kg^2.m^-1.s^-2
!               = --------------------
!                      kg.m^2.s^-2
!
!               = 0.1 kg.m^-3
!
!       Thus the result is scaled by 0.1 to return density in units
!       of kg.m^-3.
!
! COMMENTS:
!       Note that the rank-1 interface of this routine can be called with
!       either a scalar molecular weight (e.g. that for dry air for a number
!       of pressure levels) or a rank-1 molecular weight (e.g. that for air
!       where the water vapour contribution has been taken into account.)
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 01-May-2000
!                       paul.vandelst@ssec.wisc.edu
!S-
!--------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION Density_scalar( Pressure,         &  ! Input
                           Temperature,      &  ! Input
                           Molecular_Weight, &  ! Input
                           Message_Log )     &  ! Error messaging
                         RESULT( Rho )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Temperature
    REAL(fp),               INTENT(IN) :: Molecular_Weight
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Rho
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Density'
    REAL(fp),     PARAMETER :: SCALE_FACTOR = 0.1_fp
    
    ! Setup
    Rho = -ONE
    IF ( Pressure         < TOLERANCE .OR. &
         Temperature      < TOLERANCE .OR. &
         Molecular_Weight < TOLERANCE ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input Pressures/Temperature/MW < or = 0.0 found.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate density
    Rho = SCALE_FACTOR * Pressure * Molecular_Weight / ( R0 * Temperature ) 

  END FUNCTION Density_scalar


  ! ===========================================
  ! Rank-1 version with scalar molecular weight
  ! ===========================================
  FUNCTION Density_rank1a( Pressure,         &  ! Input
                           Temperature,      &  ! Input
                           Molecular_Weight, &  ! Input
                           Message_Log )     &  ! Error messaging
                         RESULT( Rho )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),               INTENT(IN) :: Molecular_Weight
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Rho(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Density'
    ! Local variables
    INTEGER :: i, n
    
    ! Setup
    Rho = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Pressure/Temperature array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate density
    DO i = 1, n
      Rho( i ) = Density_scalar( Pressure( i ), &
                                 Temperature( i ), &
                                 Molecular_Weight, &
                                 Message_Log=Message_Log )
      IF ( Rho( i ) < ZERO ) RETURN
    END DO

  END FUNCTION Density_rank1a


  ! ===========================================
  ! Rank-1 version with rank-1 molecular weight
  ! ===========================================
  FUNCTION Density_rank1b( Pressure,         &  ! Input
                           Temperature,      &  ! Input
                           Molecular_Weight, &  ! Input
                           Message_Log )     &  ! Error messaging
                         RESULT( Rho )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),               INTENT(IN) :: Molecular_Weight(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Rho(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'DENSITY'
    ! Local variables
    INTEGER :: i, n
    
    ! Setup
    Rho = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature)      /= n .OR. &
         SIZE(Molecular_Weight) /= n      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate density
    DO i = 1, n
      Rho( i ) = Density_scalar( Pressure( i ), &
                                 Temperature( i ), &
                                 Molecular_Weight( i ), &
                                 Message_Log=Message_Log )
      IF ( Rho( i ) < ZERO ) RETURN
    END DO

  END FUNCTION Density_rank1b


!--------------------------------------------------------------------------------
!S+
! NAME:
!       SVP_Water
!
! PURPOSE:
!       Function to calculate the saturation vapor pressure over water.
!
! CALLING SEQUENCE:
!       svp = SVP_Water( Temperature,            &  ! Input
!                        Message_Log=Message_Log )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Temperature:      Temperatures for which the saturation vapor
!                         pressure is required.
!                         UNITS:      Kelvin, K
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar or Rank-1
!                         ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:      Character string specifying a filename in which any
!                         Messages will be logged. If not specified, or if an
!                         error occurs opening the log file, the default action
!                         is to output Messages to standard output.
!                         UNITS:      N/A
!                         TYPE:       CHARACTER(*)
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! FUNCTION RESULT:
!       svp:               The saturation vapor pressure over water
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Temperature
!
! RESTRICTIONS:
!       Valid temperature range is 188K - 343K (-85C - +70C). A warning is
!       reported if the input temperatures are outside this range.
!
! PROCEDURE:
!       Flatau,P.J., R.L.Walko, and W.R.Cotton, 1992: "Polynomial fits to
!         saturation vapor pressure", J.Appl.Met., v31, pp1507-1513
!
!                           __ N
!                          \            i
!         SVP_Water = c0 +  >   c(i) . T
!                          /__ 
!                             i=1
!
!       where the c(i) are the relative error norm coefficients obtained
!       from the reference above.
!
!       Horner's method is used to evaluate the above polynomial.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 24-Apr-1998
!                       paul.vandelst@ssec.wisc.edu
!S-
!--------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION SVPw_scalar( Temperature,   &  ! Input
                        Message_Log  ) &  ! Error messaging
                      RESULT( svp )
    ! Arguments
    REAL(fp),          INTENT(IN) :: Temperature
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: svp
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'SVP_Water'
    INTEGER,  PARAMETER :: N_COEFFICIENTS = 8
    REAL(fp), PARAMETER :: COEFFICIENTS(0:N_COEFFICIENTS) = &
      (/ 6.11583699e+00_fp, 4.44606896e-01_fp, 1.43177157e-02_fp, &
         2.64224321e-04_fp, 2.99291081e-06_fp, 2.03154182e-08_fp, &
         7.02620698e-11_fp, 3.79534310e-14_fp,-3.21582393e-16_fp /)
    REAL(fp), PARAMETER :: MIN_TEMPERATURE = 188.15_fp
    REAL(fp), PARAMETER :: MAX_TEMPERATURE = 343.15_fp
    ! Local variables
    CHARACTER(256) :: Message
    INTEGER :: i
    REAL(fp) :: T_Celsius

    ! Setup
    IF ( Temperature < MIN_TEMPERATURE .OR. &
         Temperature > MAX_TEMPERATURE ) THEN
      WRITE( Message,'("Input Temperature ",f6.2,&
                      &" outside valid range: ",&
                      &f6.2,"K < T < ",f6.2,"K" )' ) &
                      Temperature, MIN_TEMPERATURE, MAX_TEMPERATURE
      CALL Display_Message( ROUTINE_NAME, &
                            TRIM(Message), &
                            WARNING, &
                            Message_Log=Message_Log )
    ENDIF

    ! Calculate saturation vapor pressure
    T_Celsius = Temperature - CELSIUS_TO_KELVIN
    svp       = COEFFICIENTS(N_COEFFICIENTS)
    DO i = N_COEFFICIENTS-1, 0, -1
      svp = ( svp * T_Celsius ) + COEFFICIENTS(i)
    END DO
 
  END FUNCTION SVPw_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION SVPw_rank1( Temperature,   &  ! Input
                       Message_Log  ) &  ! Error messaging
                     RESULT( svp )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Temperature(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: svp(SIZE(Temperature))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'SVP_Water'
    ! Local variables
    INTEGER :: i

    ! Calculate saturation vapor pressure
    DO i = 1, SIZE(Temperature)
      SVP( i ) = SVPw_scalar( Temperature( i ), &
                              Message_Log=Message_Log )
    END DO

  END FUNCTION SVPw_rank1


!--------------------------------------------------------------------------------
!S+
! NAME:
!       SVP_Ice
!
! PURPOSE:
!       Function to calculate the saturation vapor pressure over ice
!
! CALLING SEQUENCE:
!       svp = SVP_Ice( Temperature,            &  ! Input
!                      Message_Log=Message_Log )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Temperature:      Temperatures for which the saturation vapor
!                         pressure is required.
!                         UNITS:      Kelvin, K
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar or Rank-1
!                         ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:      Character string specifying a filename in which any
!                         Messages will be logged. If not specified, or if an
!                         error occurs opening the log file, the default action
!                         is to output Messages to standard output.
!                         UNITS:      N/A
!                         TYPE:       CHARACTER(*)
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! FUNCTION RESULT:
!       svp:               The saturation vapor pressure over ice
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Temperature
!
! RESTRICTIONS:
!       Valid Temperature range is 183K - 273K (-90C - 0C). An warning is
!       reported if the input Temperatures are outside this range.
!
! PROCEDURE:
!       Flatau,P.J., R.L.Walko, and W.R.Cotton, 1992: "Polynomial fits to
!         saturation vapor pressure", J.Appl.Met., v31, pp1507-1513
!
!                         __ N
!                        \            i
!         SVP_Ice = c0 +  >   c(i) . T
!                        /__ 
!                           i=1
!
!       where the c(i) are the relative error norm coefficients obtained
!       from the reference above.
!
!       Horner's method is used to evaluate the above polynomial.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 24-Apr-1998
!                       paul.vandelst@ssec.wisc.edu
!S-
!--------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION SVPi_scalar( Temperature,   &  ! Input
                        Message_Log  ) &  ! Error messaging
                      RESULT( svp )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Temperature
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: svp
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'SVP_Ice'
    INTEGER,  PARAMETER :: N_COEFFICIENTS = 8
    REAL(fp), PARAMETER :: COEFFICIENTS(0:N_COEFFICIENTS) = &
      (/ 6.09868993e+00_fp, 4.99320233e-01_fp, 1.84672631e-02_fp, &
         4.02737184e-04_fp, 5.65392987e-06_fp, 5.21693933e-08_fp, &
         3.07839583e-10_fp, 1.05785160e-12_fp, 1.61444444e-15_fp /)
    REAL(fp), PARAMETER :: MIN_TEMPERATURE = 183.15_fp
    REAL(fp), PARAMETER :: MAX_TEMPERATURE = 273.15_fp
    ! Local variables
    CHARACTER(256) :: Message
    INTEGER  :: i
    REAL(fp) :: T_Celsius
    
    ! Setup
    IF ( Temperature < MIN_TEMPERATURE .OR. &
         Temperature > MAX_TEMPERATURE ) THEN
      WRITE( Message,'("Input Temperature ",f6.2,&
                      &" outside valid range: ",&
                      &f6.2,"K < T < ",f6.2,"K")' ) &
                      Temperature, MIN_TEMPERATURE, MAX_TEMPERATURE
      CALL Display_Message( ROUTINE_NAME, &
                            TRIM(Message), &
                            WARNING, &
                            Message_Log=Message_Log )
    ENDIF

    ! Calculate saturation vapor pressure
    T_Celsius = Temperature - CELSIUS_TO_KELVIN
    svp       = COEFFICIENTS(N_COEFFICIENTS)
    DO i = N_COEFFICIENTS-1, 0, -1
      svp = ( svp * T_Celsius ) + COEFFICIENTS(i)
    END DO
 
  END FUNCTION SVPi_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION SVPi_rank1( Temperature,   &  ! Input
                       Message_Log  ) &  ! Error messaging
                     RESULT( svp )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Temperature(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: svp(SIZE(Temperature))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'SVP_Ice'
    ! Local variables
    INTEGER :: i

    ! Calculate saturation vapor pressure
    DO i = 1, SIZE(Temperature)
      SVP( i ) = SVPi_scalar( Temperature( i ), &
                              Message_Log=Message_Log )
    END DO

  END FUNCTION SVPi_rank1


!--------------------------------------------------------------------------------
!S+
! NAME:
!       Saturation_Mixing_Ratio
!
! PURPOSE:
!       Function to calculate the saturation mixing ratio for a given
!       pressure and temperature
!
! CALLING SEQUENCE:
!       smr = Saturation_Mixing_Ratio( Pressure,                        &  ! Input
!                                      Temperature,                     &  ! Input
!                                      Ice_Temperature=Ice_Temperature, &  ! optional input
!                                      Min_Pressure   =Min_Pressure,    &  ! Optional input
!                                      Message_Log    =Message_Log      )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Total atmospheric pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1
!                          ATTRIBUTES: INTENT(IN)
!
!       Temperature:       Atmospheric Temperature.
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Ice_Temperature:   Temperature below which the saturation vapor
!                          pressure over ice is used in the conversion.
!                          By default, only the saturation vapor pressure
!                          over water is used.
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!   
!       Min_Pressure:      Pressure value below which the saturation
!                          mixing ratio is not calculated. The default
!                          is 50hPa. Saturation mixing ratios at pressures
!                          less than the minimum pressure are set to zero.
!                          This is because at pressures less than 50mb, the
!                          saturation vapour pressure, which is based only on
!                          temperature, can exceed the total air pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!   
!       Message_Log:       Character string specifying a filename in which any
!                          Messages will be logged. If not specified, or if an
!                          error occurs opening the log file, the default action
!                          is to output Messages to standard output.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! FUNCTION RESULT:
!       smr:               The saturation mixing ratio
!                          If an error occurs, -1.0 is returned.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       The saturation mixing ratio can be defined as:
!
!                rho_ws
!          ws = --------     .....(1)
!                rho_d
!
!       where rho_ws = the partial density of water vapour required to 
!                      saturate air with respect to water at a Temperature, T
!             rho_d  = the partial density of dry air.
!
!       Equation (1) can be rewritten as:
!
!                   es
!               ---------
!                R_w . T
!         ws = ------------
!                p - es
!               ---------
!                R_d . T
!
!               R_d       es
!            = ----- . --------
!               R_w     p - es
!
!               M_w       es
!            = ----- . --------     .....(2)
!               M_d     p - es
!
!       where M_w = molecular weight of water
!             M_d = molecular weight of dry air
!             es  = water vapor partial pressure
!             p   = total air pressure
!             R_d = gas constant for dry air
!             R_w = gas constant for water vapor
!
!       The units of equation (2) are:
!
!               g     hPa
!         ws = --- . -----
!               g     hPa
!
!                      g
!            = 1000.0 ----
!                      kg
!
!       A factor of 1000 is used to return values in units of g/kg.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 18-Mar-1998
!                       paul.vandelst@ssec.wisc.edu
!S-
!--------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION SMR_scalar( Pressure,        &  ! Input
                       Temperature,     &  ! Input
                       Ice_Temperature, &  ! Optional Input
                       Min_Pressure,    &  ! Optional Input
                       Message_Log )    &  ! Error messaging
                     RESULT( smr )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Temperature
    REAL(fp),     OPTIONAL, INTENT(IN) :: Ice_Temperature
    REAL(fp),     OPTIONAL, INTENT(IN) :: Min_Pressure
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: smr
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Saturation_Mixing_Ratio'
    REAL(fp),     PARAMETER :: DEFAULT_MIN_PRESSURE = 50.0_fp
    ! Local variables
    LOGICAL :: Use_Ice_Temperature
    REAL(fp) :: Ice_T
    REAL(fp) :: Min_P
    REAL(fp) :: svp
    REAL(fp) :: dp

    ! Setup
    smr = -ONE
    IF ( Pressure    < TOLERANCE .OR. &
         Temperature < TOLERANCE ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input argument values < or = 0.0 found.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( PRESENT( Min_Pressure ) ) THEN
      Min_P = Min_Pressure
    ELSE
      Min_P = DEFAULT_MIN_PRESSURE
    END IF
    ! Default is NOT to use ice Temperature
    Use_Ice_Temperature = .FALSE.
    Ice_T               = ZERO
    ! ...unless Ice_Temperature argument is present
    IF ( PRESENT( Ice_Temperature ) ) THEN
      Use_Ice_Temperature = .TRUE.
      Ice_T               = Ice_Temperature
    END IF

    ! Only do calculation at pressures greater than the minumum
    Pressure_check: IF ( Pressure > Min_P ) THEN

      ! Calculate saturation vapor pressure
      IF ( Use_Ice_Temperature .AND. ( Temperature < Ice_T ) ) THEN
        svp = SVPi_scalar( Temperature, Message_Log=Message_Log )
      ELSE
        svp = SVPw_scalar( Temperature, Message_Log=Message_Log )
      END IF

      ! Calculate saturation mixing ratio only if the
      ! total pressure is greater than the saturation
      ! vapor pressure.
      dp = Pressure - svp
      IF ( dp > ZERO ) THEN
        smr = KG_TO_G * EPS * svp / dp
      ELSE
        smr = ZERO
      END IF

    ELSE

      smr = ZERO

    END IF Pressure_check

  END FUNCTION SMR_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION SMR_rank1( Pressure,        &  ! Input
                      Temperature,     &  ! Input
                      Ice_Temperature, &  ! Optional Input
                      Min_Pressure,    &  ! Optional Input
                      Message_Log )    &  ! Error messaging
                    RESULT( smr )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Ice_Temperature
    REAL(fp),     OPTIONAL, INTENT(IN) :: Min_Pressure
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: smr(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Saturation_Mixing_Ratio'
    ! Local variables
    INTEGER :: i, n

    ! Setup
    smr = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature) /=  n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Pressure/Temperature array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate saturation mixing ratio
    DO i = 1, n
      smr( i ) = SMR_scalar( Pressure( i ), &
                             Temperature( i ), &
                             Ice_Temperature = Ice_Temperature, &
                             Min_Pressure    = Min_Pressure, &
                             Message_Log     = Message_Log )
      IF ( smr( i ) < ZERO ) RETURN
    END DO

  END FUNCTION SMR_rank1


!--------------------------------------------------------------------------------
!S+
! NAME:
!       Virtual_Temperature
!
! PURPOSE:
!       Function to calculate the virtual temperature given the temperature and
!       water vapor mixing ratio.
!
!       There is also an "inverse" capability to compute the temperature given
!       the virtual temperature and the water vapor mixing ratio.
!
! CALLING SEQUENCE:
!       Tv = Virtual_Temperature( Temperature,              &  ! Input
!                                 Water_Vapor_Mixing_Ratio, &  ! Input
!                                 Inverse    =Inverse,      &  ! Optional Input
!                                 Message_Log=Message_Log   )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Temperature:               Atmospheric temperature.
!                                  UNITS:      Kelvin, K
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar or Rank-1
!                                  ATTRIBUTES: INTENT(IN)
!   
!       Water_Vapor_Mixing_Ratio:  Water vapor mass mixing ratio.
!                                  UNITS:      g/kg
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Same as Temperature
!                                  ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Inverse:                   Set this argument to compute the temperature
!                                  given the input is the virtual temperature.
!                                  If = 0, virtual temperature is calculated (default)
!                                     = 1, input value is assumed to be the virtual
!                                          temperature and the temperature is
!                                          calculated.
!                                  UNITS:      N/A
!                                  TYPE:       INTEGER
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: OPTIONAL, INTENT(IN)
!
!       Message_Log:               Character string specifying a filename in which any
!                                  Messages will be logged. If not specified, or if an
!                                  error occurs opening the log file, the default action
!                                  is to output Messages to standard output.
!                                  UNITS:      N/A
!                                  TYPE:       CHARACTER(*)
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! FUNCTION RESULT:
!       Tv:                        The virtual temperature. If the optional Inverse
!                                  argument is set, the function result is the
!                                  atmospheric temperature.
!                                  If an error occurs, -1.0 is returned.
!                                  UNITS:      Kelvin, K
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Same as input Temperature
!
! PROCEDURE:
!       The virtual temperature, the temperature that dry air must have in
!       order to have the same density as moist air at the same pressure, is
!       calculated using:
!
!                  [    eps + w    ]
!         Tv = T * [ ------------- ]     .......................(1)
!                  [ eps ( 1 + w ) ]
!
!       where T   = temperature,
!             w   = water vapour mixing ratio, and
!             eps = ratio of the molecular weights of water and dry air.
!
!       An approximation to eqn.(1) is,
!
!                  [      1 - eps    ]
!         Tv = T * [ 1 + --------- w ]
!                  [        eps      ]
!
!            = T * [ 1 + ( 0.608 * w ) ]   .....................(2)
!
!       however, depending on what accuracy is required (keeping in mind that
!       water vapor measurements are probably good to 2-5%), eqn.(2) can 
!       differ from (1) by around 0.06-0.08K near the surface.
!
!       If virtual temperature is used to calculate geopotential heights,
!       this difference can lead to errors of up to 0.6-0.7m.
!
!       So I took the slightly more computationally expensive road
!       and use eqn.(1).
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 01-May-2000
!                       paul.vandelst@ssec.wisc.edu
!S-
!--------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION Tv_scalar( Temperature,              &  ! Input
                      Water_Vapor_Mixing_Ratio, &  ! Input
                      Inverse,                  &  ! Optional input
                      Message_Log )             &  ! Error messaging
                    RESULT( Tv )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Temperature
    REAL(fp),               INTENT(IN) :: Water_Vapor_Mixing_Ratio
    INTEGER,      OPTIONAL, INTENT(IN) :: Inverse
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Tv
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Virtual_Temperature'
    ! Local variables
    LOGICAL :: Forward
    REAL(fp) :: mr_H2O
    REAL(fp) :: Wterm
    
    ! Setup
    Tv = -ONE
    IF ( Temperature < TOLERANCE .OR. Water_Vapor_Mixing_Ratio < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input temperature/water vapor mixing ratio < or = 0.0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    ! Always do the forward calculation....
    Forward = .TRUE.
    ! ...unless the Inverse argument is set
    IF ( PRESENT( Inverse ) ) THEN
      IF( Inverse == SET ) Forward = .FALSE.
    END IF

    ! Calculate the virtual temperature
    mr_H2O = G_TO_KG * Water_Vapor_Mixing_Ratio
    Wterm  = ( EPS + mr_H2O ) / ( EPS * ( ONE + mr_H2O ) )
    IF ( Forward ) THEN
      Tv = Temperature * Wterm
    ELSE
      Tv = Temperature / Wterm
    END IF

  END FUNCTION Tv_scalar 


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION Tv_rank1( Temperature,              &  ! Input
                     Water_Vapor_Mixing_Ratio, &  ! Input
                     Inverse,                  &  ! Optional input
                     Message_Log )             &  ! Error messaging
                   RESULT( Tv )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),               INTENT(IN) :: Water_Vapor_Mixing_Ratio(:)
    INTEGER,      OPTIONAL, INTENT(IN) :: Inverse
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Tv(SIZE(Temperature))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Virtual_Temperature'
    ! Local variables
    INTEGER :: i, n
    
    ! Setup
    Tv = -ONE
    n = SIZE(Temperature)
    IF ( SIZE(Water_Vapor_Mixing_Ratio) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Temperature/Water_Vapor_'//&
                            'Mixing_Ratio array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Compute the virtual temperature (or inverse)
    DO i = 1, n
      Tv( i ) = Tv_scalar( Temperature( i ), &
                           Water_Vapor_Mixing_Ratio( i ), &
                           Inverse    =Inverse, &
                           Message_Log=Message_Log )
      IF ( Tv( i ) < ZERO ) RETURN
    END DO

  END FUNCTION Tv_rank1


!--------------------------------------------------------------------------------
!S+
! NAME:
!       Potential_Temperature
!
! PURPOSE:
!       Function to calculate the potential temperature given the temperature
!       and pressure.
!
!       There is also an "inverse" capability to compute the temperature given
!       the potential temperature and pressure.
!
! CALLING SEQUENCE:
!       Theta = Potential_Temperature( Temperature,            &  ! Input
!                                      Pressure,               &  ! Input
!                                      Inverse    =Inverse,    &  ! Optional input
!                                      Message_Log=Message_Log )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Temperature:       Atmospheric temperature.
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1
!                          ATTRIBUTES: INTENT(IN)
!   
!       Pressure:          Atmospheric pressure
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Temperature
!                          ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Inverse:           Set this argument to compute the temperature
!                          given the input is the potential temperature.
!                          If = 0, potential temperature is calculated (default)
!                             = 1, input value is assumed to be the potential
!                                  temperature and the temperature is
!                                  calculated.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!
!       Message_Log:       Character string specifying a filename in which any
!                          Messages will be logged. If not specified, or if an
!                          error occurs opening the log file, the default action
!                          is to output Messages to standard output.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! FUNCTION RESULT:
!       Theta:             The potential temperature. If the optional Inverse
!                          argument is set, the function result is the
!                          atmospheric temperature.
!                          If an error occurs, -1.0 is returned.
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Temperature
!
! PROCEDURE:
!       The potential temperature of a parcel of air is that temperature
!       the parcel would have if it were expanded or compressed adiabatically
!       to some reference pressure.
!
!       The conversion is given by Poisson's equation:
!
!                             R/Cp
!                     [  P0  ]
!         Theta = T * [ ---- ]           .......................(1)
!                     [  p   ]
!
!       where T   = temperature,
!             p   = pressure,
!             P0  = standard pressure
!             R   = gas constant,
!             Cp  = specific heat of gas at constant pressure.
!
!       This routine uses the standard atmosphere as the reference pressure,
!       R for dry air, and Cp for an ideal diatomic gas:
!
!               7      R0
!         Cp = --- . -------
!               2     MWair
!
!       where R0    = universal gas constant
!             MWair = molecular weight of dry air
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 11-Sep-2002
!                       paul.vandelst@ssec.wisc.edu
!S-
!--------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION Theta_scalar( Temperature,  &  ! Input
                         Pressure,     &  ! Input
                         Inverse,      &  ! Optional Input
                         Message_Log ) &  ! Error messaging
                       RESULT( Theta )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Temperature
    REAL(fp),               INTENT(IN) :: Pressure
    INTEGER,      OPTIONAL, INTENT(IN) :: Inverse
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Theta
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Potential_Temperature'
    REAL(fp),     PARAMETER :: EXPONENT_TERM = R_DRYAIR / Cp_DRYAIR
    ! Local variables
    LOGICAL  :: Forward
    REAL(fp) :: PressureTerm
    
    ! Setup
    Theta = -ONE
    IF ( Temperature < TOLERANCE .OR. Pressure < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input temperature/pressure < or = 0.0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    ! Always do the forward calculation....
    Forward = .TRUE.
    ! ...unless the Inverse argument is set
    IF ( PRESENT( Inverse ) ) THEN
      IF( Inverse == SET ) Forward = .FALSE.
    END IF

    ! Calculate the potential temperature
    PressureTerm = ( ( P0 / Pressure )**EXPONENT_TERM )
    IF ( Forward ) THEN
      Theta = Temperature * PressureTerm
    ELSE
      Theta = Temperature / PressureTerm
    END IF

  END FUNCTION Theta_scalar 


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION Theta_rank1( Temperature,  &  ! Input
                        Pressure,     &  ! Input
                        Inverse,      &  ! Optional Input
                        Message_Log ) &  ! Error messaging
                      RESULT( Theta )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),               INTENT(IN) :: Pressure(:)
    INTEGER,      OPTIONAL, INTENT(IN) :: Inverse
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Theta(SIZE(Temperature))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Potential_Temperature'
    ! Local variables
    INTEGER :: i, n

    ! Setup
    Theta = -ONE
    n = SIZE(Temperature)
    IF ( SIZE(Pressure) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Temperature/Pressure array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate the potential temperature
    DO i = 1, n
      Theta( i ) = Theta_scalar( Temperature( i ), &
                                 Pressure( i ), &
                                 Inverse    =Inverse, &
                                 Message_Log=Message_Log )
      IF ( Theta( i ) < ZERO ) RETURN
    END DO

  END FUNCTION Theta_rank1

END MODULE Atmospheric_Properties
