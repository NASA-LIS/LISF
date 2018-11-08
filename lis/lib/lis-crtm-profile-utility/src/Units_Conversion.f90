!
! Units_Conversion
!
! Module containing routines to convert atmospheric profile
! concentration units.
!
!
! CONTAINS:
!       Specific amount <-> mixing ratio
!       --------------------------------
!         SA_to_MR:       Function to convert gas concentrations from
!                         specific amount to mixing ratio.
!
!         MR_to_SA:       Function to convert gas concentrations from
!                         mixing ratio to specific amount.
!
!
!       Relative humidity <-> mixing ratio
!       ----------------------------------
!         RH_to_MR:       Function to convert water vapor amounts from
!                         relative humidity to mixing ratio.
!
!         MR_to_RH:       Function to convert water vapor amounts from
!                         mixing ratio to relative humidity.
!
!
!       Mixing ratio <-> ppmv
!       ---------------------
!         MR_to_PPMV:     Function to convert gas concentrations from 
!                         mixing ratio to parts-per-million by volume.
!
!         PPMV_to_MR:     Function to convert gas concentrations from
!                         parts-per-million by volume to mixing ratio.
!
!
!       Partial pressure <-> ppmv
!       -------------------------
!         PPMV_to_PP:     Function to convert gas concentrations from
!                         parts-per-million by volume to partial pressure.
!
!         PP_to_PPMV:     Function to convert gas concentrations from
!                         partial pressure to parts-per-million by volume.
!
!
!       Partial pressure <-> mixing ratio
!       ---------------------------------
!         MR_to_PP:       Function to convert gas concentrations from mixing
!                         ratio to partial pressure in hectoPascals
!
!         PP_to_MR:       Function to convert gas concentrations from partial
!                         pressure in hectoPascals to mixing ratio in g/kg.
!
!
!       Partial pressure <-> mass density
!       ---------------------------------
!         PP_to_MD:       Function to convert gas concentrations in pressure
!                         units to mass density.
!
!         MD_to_PP:       Function to convert gas concentration mass density
!                         to partial pressure.
!
!
!       Partial pressure <-> number density
!       -----------------------------------
!         PP_to_ND:       Function to convert gas concentrations from (partial)
!                         pressures in hectoPascals to molecules/m^3.
!
!         ND_to_PP:       Function to convert gas concentrations from number
!                         densities in molecules/m^3 to (partial) pressures
!                         in hectoPascals.
!
!
!       Number density <-> ppmv
!       -----------------------
!         PPMV_to_ND:     Function to convert gas concentrations from
!                         parts-per-million by volume to molecules/m^3.
!
!         ND_to_PPMV:     Function to convert gas concentrations from number
!                         densities in molecules/m^3 to parts-per-million
!                         by volume.
!
!
!       Column density <-> ppmv
!       -----------------------
!         PPMV_to_KMOL:   Function to convert gas concentrations from parts-per-
!                         million by volume to kilomoles per cm^2.
!
!         KMOL_to_PPMV:   Function to convert gas concentrations from kilomoles
!                         per cm^2 to parts-per-million by volume.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 01-May-2000
!                       paul.vandelst@noaa.gov
!

MODULE Units_Conversion

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Modules used
  USE Type_Kinds, ONLY: fp
  USE Message_Handler
  USE Atmospheric_Properties
  USE Profile_Utility_Parameters
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: SA_to_MR,     MR_to_SA
  PUBLIC :: RH_to_MR,     MR_to_RH
  PUBLIC :: MR_to_PPMV,   PPMV_to_MR
  PUBLIC :: PPMV_to_PP,   PP_to_PPMV
  PUBLIC :: MR_to_PP,     PP_to_MR
  PUBLIC :: PP_to_MD,     MD_to_PP
  PUBLIC :: PP_to_ND,     ND_to_PP
  PUBLIC :: PPMV_to_ND,   ND_to_PPMV
  PUBLIC :: PPMV_to_KMOL, KMOL_to_PPMV


  ! ---------------------
  ! Procedure overloading
  ! ---------------------
  INTERFACE SA_to_MR
    MODULE PROCEDURE SA_to_MR_scalar
    MODULE PROCEDURE SA_to_MR_rank1
  END INTERFACE SA_to_MR

  INTERFACE MR_to_SA
    MODULE PROCEDURE MR_to_SA_scalar
    MODULE PROCEDURE MR_to_SA_rank1
  END INTERFACE MR_to_SA

  INTERFACE RH_to_MR
    MODULE PROCEDURE RH_to_MR_scalar
    MODULE PROCEDURE RH_to_MR_rank1
  END INTERFACE RH_to_MR

  INTERFACE MR_to_RH
    MODULE PROCEDURE MR_to_RH_scalar
    MODULE PROCEDURE MR_to_RH_rank1
  END INTERFACE MR_to_RH

  INTERFACE MR_to_PPMV
    MODULE PROCEDURE MR_to_PPMV_scalar
    MODULE PROCEDURE MR_to_PPMV_rank1
  END INTERFACE MR_to_PPMV

  INTERFACE PPMV_to_MR
    MODULE PROCEDURE PPMV_to_MR_scalar
    MODULE PROCEDURE PPMV_to_MR_rank1
  END INTERFACE PPMV_to_MR

  INTERFACE PPMV_to_PP
    MODULE PROCEDURE PPMV_to_PP_scalar
    MODULE PROCEDURE PPMV_to_PP_rank1
  END INTERFACE PPMV_to_PP

  INTERFACE PP_to_PPMV
    MODULE PROCEDURE PP_to_PPMV_scalar
    MODULE PROCEDURE PP_to_PPMV_rank1
  END INTERFACE PP_to_PPMV

  INTERFACE MR_to_PP
    MODULE PROCEDURE MR_to_PP_scalar
    MODULE PROCEDURE MR_to_PP_rank1
  END INTERFACE MR_to_PP

  INTERFACE PP_to_MR
    MODULE PROCEDURE PP_to_MR_scalar
    MODULE PROCEDURE PP_to_MR_rank1
  END INTERFACE PP_to_MR

  INTERFACE PP_to_MD
    MODULE PROCEDURE PP_to_MD_scalar
    MODULE PROCEDURE PP_to_MD_rank1
  END INTERFACE PP_to_MD

  INTERFACE MD_to_PP
    MODULE PROCEDURE MD_to_PP_scalar
    MODULE PROCEDURE MD_to_PP_rank1
  END INTERFACE MD_to_PP

  INTERFACE PP_to_ND
    MODULE PROCEDURE PP_to_ND_scalar
    MODULE PROCEDURE PP_to_ND_rank1
  END INTERFACE PP_to_ND

  INTERFACE ND_to_PP
    MODULE PROCEDURE ND_to_PP_scalar
    MODULE PROCEDURE ND_to_PP_rank1
  END INTERFACE ND_to_PP

  INTERFACE PPMV_to_ND
    MODULE PROCEDURE PPMV_to_ND_scalar
    MODULE PROCEDURE PPMV_to_ND_rank1
  END INTERFACE PPMV_to_ND

  INTERFACE ND_to_PPMV
    MODULE PROCEDURE ND_to_PPMV_scalar
    MODULE PROCEDURE ND_to_PPMV_rank1
  END INTERFACE ND_to_PPMV

  INTERFACE PPMV_to_KMOL
    MODULE PROCEDURE PPMV_to_KMOL_scalar
    MODULE PROCEDURE PPMV_to_KMOL_rank1
  END INTERFACE PPMV_to_KMOL

  INTERFACE KMOL_to_PPMV
    MODULE PROCEDURE KMOL_to_PPMV_scalar
    MODULE PROCEDURE KMOL_to_PPMV_rank1
  END INTERFACE KMOL_to_PPMV


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_RCS_ID = &
  '$Id: Units_Conversion.f90 2018 2008-04-01 20:12:33Z paul.vandelst@noaa.gov $'
  REAL(fp), PARAMETER :: HUNDRED      = 100.0_fp
  REAL(fp), PARAMETER :: FROM_PERCENT = ONE / HUNDRED


CONTAINS


!------------------------------------------------------------------------------
!S+
! NAME:
!       SA_to_MR
!
! PURPOSE:
!       Function to convert gas specific amounts to mass mixing ratio
!
! CALLING SEQUENCE:
!       Mixing_Ratio = SA_to_MR( Specific_Amount,         &  ! Input
!                                Water_Vapor=Water_Vapor, &  ! Optional input
!                                Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Specific_Amount:   Gas specific amount comcentration.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Water_Vapor:       Water vapor mass mixing ratio. If this argument is
!                          not supplied, the mandatory input SPECIFIC_AMOUNT
!                          argument is assumed to be water vapor.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
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
!       Mixing_Ratio:      The gas mass mixing ratio. If an error occurs,
!                          the value -1.0 is returned.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Specific_Amount
!
! PROCEDURE:
!       The specific amount is defined as the ratio of the mass of gas
!       to the total mass of air and water vapor in a given volume:
!
!                  Mg
!         SA = ---------   .......(1)
!               Md + Mw
!
!       where Mg = mass of gas
!             Mw = mass of water vapor
!             Md = mass of dry air (including Mg).
!
!       The gas mass mixing ratio is defined as the ratio of the mass
!       of the gas to the mass of dry air in a given volume:
!
!                Mg
!         MRg = ----    .......(2)
!                Md
!
!       Rearranging (1) and substituting in (2) gives,
!
!          1      Md + Mw
!         ---- = ---------
!          SA       Mg
!
!                 Md     Mw
!              = ---- + ----
!                 Mg     Mg
!
!                 Md     Mw   Md                                      Md
!              = ---- + ----.----    ( multiplying the second term by ---- )
!                 Mg     Md   Mg                                      Md
!
!                  1 
!              = ----- ( 1 + MRw )
!                 MRg
!
!       therefore,
!
!         MRg = SA ( 1 + MRw )
!
!       for input units of g/g or kg/kg. For input units of g/kg then,
!
!         MRg = SA ( 1 + 0.001*MRw )   .......(3)
!
!       If the input specific amount is for water vapor (specific humidity) 
!       then (3) becomes,
!
!         MRw = SA ( 1 + 0.001*MRw )
!
!       i.e.
!
!         MRw = SA + 0.001*MRw*SA
!
!       i.e.
!
!         MRw( 1 - 0.001*SA ) = SA
!
!       therefore,
!
!                      SA
!         MRw = ------------------   .......(4)
!                ( 1 - 0.001*SA )
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 05-Dec-2003
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION SA_to_MR_scalar( Specific_Amount, &  ! Input
                            Water_Vapor,     &  ! Optional Input
                            Message_Log )    &  ! Error messaging
                          RESULT( Mixing_Ratio )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Specific_Amount
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Mixing_Ratio
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'SA_to_MR'

    ! Setup
    Mixing_Ratio = -ONE
    IF ( Specific_Amount < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input specific amount < 0.0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate mixing ratio in g/kg
    IF ( PRESENT( Water_Vapor ) ) THEN
      IF ( Water_Vapor < ZERO ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Input Water_Vapor mixing ratio < 0.0.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Mixing_Ratio = Specific_Amount * ( ONE + ( G_TO_KG * Water_Vapor ) )
    ELSE
      Mixing_Ratio = Specific_Amount / ( ONE - ( G_TO_KG * Specific_Amount ) )
    END IF

  END FUNCTION SA_to_MR_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION SA_to_MR_rank1( Specific_Amount, &  ! Input
                           Water_Vapor,     &  ! Optional Input
                           Message_Log )    &  ! Error messaging
                         RESULT( Mixing_Ratio )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Specific_Amount(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Mixing_Ratio(SIZE(Specific_Amount))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'SA_to_MR'
    ! Local variables
    INTEGER :: k

    ! Setup
    Mixing_Ratio = -ONE

    ! Calculate mixing ratio in g/kg
    IF ( PRESENT(Water_Vapor) ) THEN
    
      ! Water vapour argument is present
      IF ( SIZE(Water_Vapor) /= SIZE(Specific_Amount) ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Inconsistent input Specific_Amount/Water_Vapor array sizes.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      ENDIF
      DO k = 1, SIZE(Specific_Amount)
        Mixing_Ratio( k ) = SA_to_MR_scalar( Specific_Amount( k ), &
                                             Water_Vapor=Water_Vapor(k), &
                                             Message_Log=Message_Log )
        IF ( Mixing_Ratio( k ) < ZERO ) RETURN
      END DO

    ELSE

      ! Water vapour argument is NOT present
      DO k = 1, SIZE(Specific_Amount)
        Mixing_Ratio( k ) = SA_to_MR_scalar( Specific_Amount( k ), &
                                             Message_Log=Message_Log )
        IF ( Mixing_Ratio( k ) < ZERO ) RETURN
      END DO

    END IF

  END FUNCTION SA_to_MR_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       MR_to_SA
!
! PURPOSE:
!       Function to convert gas concentrations from mass mixing ratio to
!       specific amounts.
!
! CATEGORY:
!       Profile Utility
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       Specific_Amount = MR_to_SA( Mixing_Ratio,            &  ! Input
!                                   Water_Vapor=Water_Vapor, &  ! Optional input
!                                   Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Mixing_Ratio:      Gas mass mixing ratio.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Water_Vapor:       Water vapor specific humidity. If this argument is
!                          not supplied, the mandatory input MIXING_RATIO
!                          argument is assumed to be water vapor.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
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
!       Specific_Amount:   The gas specific amount. If an error occurs,
!                          the value -1.0 is returned.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Mixing_Ratio
!
! PROCEDURE:
!       From the SA_to_MR conversion, we know that for input units of
!       g/kg that,
!
!         MRg = SAg ( 1 + 0.001*MRw )   .......(1)
!
!       and 
!                      SAw
!         MRw = -------------------   .......(2)
!                ( 1 - 0.001*SAw )
!
!       where MRg = mass mixing ratio of gas
!             MRw = mass mixing ratio of water vapor
!             SAg = specific amount of gas
!             SAw = specific amount of water vapor (specific humidity)
!
!       Rearranging (1) and (2) gives,
!
!         SAg = MRg ( 1 - 0.001*SAw )
!
!       and for water vapor only,
!
!                       MRw
!         SAw = -------------------
!                ( 1 + 0.001*MRw )
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 05-Dec-2003
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION MR_to_SA_scalar( Mixing_Ratio, &  ! Input
                            Water_Vapor,  &  ! Optional Input
                            Message_Log ) &  ! Error messaging
                          RESULT( Specific_Amount )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Mixing_Ratio
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Specific_Amount
    ! Local parameters
    CHARACTER(*),  PARAMETER :: ROUTINE_NAME = 'MR_to_SA'

    ! Setup
    Specific_Amount = -ONE
    IF ( Mixing_Ratio < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input mixing ratio < 0.0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate specific amount in g/kg
    IF ( PRESENT(Water_Vapor) ) THEN
      IF ( Water_Vapor < ZERO ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Input Water_Vapor specific humidity < 0.0.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Specific_Amount = Mixing_Ratio * ( ONE - ( G_TO_KG * Water_Vapor ) )
    ELSE
      Specific_Amount = Mixing_Ratio / ( ONE + ( G_TO_KG * Mixing_Ratio ) )
    END IF

  END FUNCTION MR_to_SA_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION MR_to_SA_rank1( Mixing_Ratio, &  ! Input
                           Water_Vapor,  &  ! Optional Input
                           Message_Log ) &  ! Error messaging
                         RESULT( Specific_Amount )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Mixing_Ratio(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Specific_Amount(SIZE(Mixing_Ratio))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MR_to_SA'
    ! Local variables
    INTEGER :: k

    ! Setup
    Specific_Amount = -ONE

    ! Calculate specific amount in g/kg
    IF ( PRESENT(Water_Vapor) ) THEN

      ! Water vapour argument is present
      IF ( SIZE(Water_Vapor) /= SIZE(Mixing_Ratio) ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Inconsistent input Mixing_Ratio/Water_Vapor array sizes.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      ENDIF
      DO k = 1, SIZE(Mixing_Ratio)

        Specific_Amount( k ) = MR_to_SA_scalar( Mixing_Ratio( k ), &
                                                Water_Vapor=Water_Vapor(k), &
                                                Message_Log=Message_Log )
        IF ( Specific_Amount( k ) < ZERO ) RETURN
      END DO

    ELSE

      ! Water vapour argument is NOT present
      DO k = 1, SIZE(Mixing_Ratio)
        Specific_Amount( k ) = MR_to_SA_scalar( Mixing_Ratio( k ), &
                                                Message_Log=Message_Log )
        IF ( Specific_Amount( k ) < ZERO ) RETURN
      END DO

    END IF

  END FUNCTION MR_to_SA_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       RH_to_MR
!
! PURPOSE:
!       Function to convert relative humidity to water vapor mixing ratio
!
! CALLING SEQUENCE:
!       Mixing_Ratio = RH_to_MR( Pressure,                        &  ! Input
!                                Temperature,                     &  ! Input
!                                Relative_Humidity,               &  ! Input
!                                Ice_Temperature=Ice_Temperature, &  ! optional input
!                                Min_Pressure   =Min_Pressure,    &  ! Optional input
!                                Message_Log    =Message_Log      )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Total atmospheric pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!
!       Temperature:       Atmospheric temperature.
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as pressure
!                          ATTRIBUTES: INTENT(IN)
!
!       Relative_Humidity: Atmospheric relative humidity.
!                          UNITS:      %
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as pressure
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
!                          is 50mb. Saturation mixing ratios below the
!                          minimum pressure are set to zero. This is
!                          because at pressures less than 50mb, the
!                          saturation vapour Pressure, which is based
!                          only on temperature, can exceed the total
!                          air Pressure.
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
!       Mixing_Ratio:      Water vapor mass mixing ratio. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       Once the saturation mixing ratio is calculated the mixing ratio
!       corresponding to the input relative humidity is determined using:
!
!                       Relative_Humidity * Saturation_Mixing_Ratio
!       Mixing_Ratio = ---------------------------------------------
!                                        100
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 02-Mar-1999
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION RH_to_MR_scalar( Pressure,          &  ! Input
                            Temperature,       &  ! Input
                            Relative_Humidity, &  ! Input
                            Ice_Temperature,   &  ! Optional Input
                            Min_Pressure,      &  ! Optional Input
                            Message_Log )      &  ! Error messaging
                          RESULT( Mixing_Ratio )
    ! Arguments
    REAL(fp),              INTENT(IN) :: Pressure
    REAL(fp),              INTENT(IN) :: Temperature
    REAL(fp),              INTENT(IN) :: Relative_Humidity
    REAL(fp),    OPTIONAL, INTENT(IN) :: Ice_Temperature
    REAL(fp),    OPTIONAL, INTENT(IN) :: Min_Pressure
    CHARACTER(*),OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Mixing_Ratio
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'RH_to_MR'
    ! Local variables
    REAL(fp) :: smr

    ! Setup
    Mixing_Ratio = -ONE
    IF ( Pressure          < TOLERANCE .OR. &
         Temperature       < TOLERANCE .OR. &
         Relative_Humidity < ZERO      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input argument values < or = 0.0 found.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( Relative_Humidity > HUNDRED ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input Relative_Humidity > 100%', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate saturation mixing ratio in g/kg
    smr = Saturation_Mixing_Ratio( Pressure, &
                                   Temperature, &
                                   Ice_Temperature=Ice_Temperature, &
                                   Min_Pressure   =Min_Pressure )
    IF ( smr < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error calculating saturation mixing ratio.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate mixing ratio in g/kg
    Mixing_Ratio = FROM_PERCENT * Relative_Humidity * smr

  END FUNCTION RH_to_MR_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION RH_to_MR_rank1( Pressure,          &  ! Input
                           Temperature,       &  ! Input
                           Relative_Humidity, &  ! Input
                           Ice_Temperature,   &  ! Optional Input
                           Min_Pressure,      &  ! Optional Input
                           Message_Log )      &  ! Error messaging
                         RESULT( Mixing_Ratio )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),               INTENT(IN) :: Relative_Humidity(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Ice_Temperature
    REAL(fp),     OPTIONAL, INTENT(IN) :: Min_Pressure
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Mixing_Ratio(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'RH_to_MR'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    Mixing_Ratio = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature)       /= n .OR. & 
         SIZE(Relative_Humidity) /= n      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate mixing ratio in g/kg
    DO k = 1, n
      Mixing_Ratio( k ) = RH_to_MR_scalar( Pressure( k ), &
                                           Temperature( k ), &
                                           Relative_Humidity( k ), &
                                           Ice_Temperature=Ice_Temperature, &
                                           Min_Pressure   =Min_Pressure, &
                                           Message_Log    =Message_Log )
      IF ( Mixing_Ratio( k ) < ZERO ) RETURN
    END DO

  END FUNCTION RH_to_MR_rank1



!------------------------------------------------------------------------------
!S+
! NAME:
!       MR_to_RH
!
! PURPOSE:
!       Function to convert water vapor mixing ratio to relative humidity
!
! CALLING SEQUENCE:
!       Relative_Humidty = MR_to_RH( Pressure,                        &  ! Input
!                                    Temperature,                     &  ! Input
!                                    Mixing_Ratio,                    &  ! Input
!                                    Ice_Temperature=Ice_Temperature, &  ! optional input
!                                    Min_Pressure   =Min_Pressure,    &  ! Optional input
!                                    Message_Log    =Message_Log      )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Total atmospheric pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!
!       Temperature:       Atmospheric temperature.
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as pressure
!                          ATTRIBUTES: INTENT(IN)
!
!       Mixing_Ratio:      Water vapor mixing ratio.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as pressure
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
!                          is 50mb. Saturation mixing ratios below the
!                          minimum pressure are set to zero. This is
!                          because at pressures less than 50mb, the
!                          saturation vapour Pressure, which is based
!                          only on temperature, can exceed the total
!                          air Pressure.
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
!       Relative_Humidity: Relative humidity. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      %
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       Once the saturation mixing ratio is calculated the relative humidity
!       corresponding to the input mixing ratio is determined using:
!
!                                         Mixing_Ratio
!       Relative_Humidity = 100.0 * -------------------------
!                                    Saturation_Mixing_Ratio
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 02-Mar-1999
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION MR_to_RH_scalar( Pressure,        &  ! Input
                            Temperature,     &  ! Input
                            Mixing_Ratio,    &  ! Input
                            Ice_Temperature, &  ! Optional Input
                            Min_Pressure,    &  ! Optional Input
                            Message_Log )    &  ! Error messaging
                          RESULT( Relative_Humidity )
    ! Arguments
    REAL(fp),           INTENT(IN) :: Pressure
    REAL(fp),           INTENT(IN) :: Temperature
    REAL(fp),           INTENT(IN) :: Mixing_Ratio
    REAL(fp), OPTIONAL, INTENT(IN) :: Ice_Temperature
    REAL(fp), OPTIONAL, INTENT(IN) :: Min_Pressure
    CHARACTER(*),  OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Relative_Humidity
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MR_to_RH'
    ! Local variables
    REAL(fp) :: smr

    ! Setup
    Relative_Humidity = -ONE
    IF ( Pressure     < TOLERANCE .OR. &
         Temperature  < TOLERANCE .OR. &
         Mixing_Ratio < ZERO           ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input argument values < or = 0.0 found.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate saturation mixing ratio in g/kg
    smr = Saturation_Mixing_Ratio( Pressure, &
                                   Temperature, &
                                   Ice_Temperature=Ice_Temperature, &
                                   Min_Pressure   =Min_Pressure, &
                                   Message_Log    =Message_Log )
    IF ( smr < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error calculating saturation mixing ratio.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate relative humidity in %
    IF ( smr > ZERO ) THEN
      Relative_Humidity = HUNDRED * Mixing_Ratio / smr
    ELSE
      Relative_Humidity = ZERO
    END IF

  END FUNCTION MR_to_RH_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION MR_to_RH_rank1( Pressure,        &  ! Input
                           Temperature,     &  ! Input
                           Mixing_Ratio,    &  ! Input
                           Ice_Temperature, &  ! Optional Input
                           Min_Pressure,    &  ! Optional Input
                           Message_Log )    &  ! Error messaging
                         RESULT( Relative_Humidity )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),               INTENT(IN) :: Mixing_Ratio(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Ice_Temperature
    REAL(fp),     OPTIONAL, INTENT(IN) :: Min_Pressure
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Relative_Humidity(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MR_to_RH'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    Relative_Humidity = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature)  /= n .OR. & 
         SIZE(Mixing_Ratio) /= n      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate relative humidity in %
    DO k = 1, n
      Relative_Humidity( k ) = MR_to_RH_scalar( Pressure( k ), &
                                                Temperature( k ), &
                                                Mixing_Ratio( k ), &
                                                Ice_Temperature=Ice_Temperature, &
                                                Min_Pressure   =Min_Pressure, &
                                                Message_Log    =Message_Log )
      IF ( Relative_Humidity( k ) < ZERO ) RETURN
    END DO

  END FUNCTION MR_to_RH_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       MR_to_PPMV
!
! PURPOSE:
!       Function to convert gas concentrations from mass mixing ratio in g/kg
!       to volume mixing ratio in ppmv
!
! CALLING SEQUENCE:
!       ppmv = MR_to_PPMV( Mixing_Ratio,            &  ! Input
!                          Molecule_ID=Molecule_ID, &  ! Input
!                          Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Mixing_Ratio:  Mass mixing ratio of gas.
!                      UNITS:      g/kg
!                      TYPE:       REAL(fp)
!                      DIMENSION:  Scalar or Rank-1 (K x 1)
!                      ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Molecule_ID:   HITRAN molecular designation identifying the
!                      molecule for which the concentration units
!                      conversion is required. If not specified, the
!                      default value is that for water vapor.
!                      Valid values are:
!                        1: H2O       9: SO2      17: HI       25: H2O2
!                        2: CO2      10: NO2      18: ClO      26: C2H2
!                        3: O3       11: NH3      19: OCS      27: C2H6
!                        4: N2O      12: HNO3     20: H2CO     28: PH3
!                        5: CO       13: OH       21: HOCl     29: COF2
!                        6: CH4      14: HF       22: N2       30: SF6
!                        7: O2       15: HCl      23: HCN      31: H2S
!                        8: NO       16: HBr      24: CH3Cl    32: HCOOH
!                      UNITS:      N/A
!                      TYPE:       INTEGER
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: OPTIONAL, INTENT(IN)
!
!       Message_Log:   Character string specifying a filename in which any
!                      Messages will be logged. If not specified, or if an
!                      error occurs opening the log file, the default action
!                      is to output Messages to standard output.
!                      UNITS:      N/A
!                      TYPE:       CHARACTER(*)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! FUNCTION RESULT:
!       ppmv:          Volume mixing ratio of gas. If an error occurs,
!                      -1.0 is returned.
!                      UNITS:      ppmv
!                      TYPE:       REAL(fp)
!                      DIMENSION:  Same as input Mixing_Ratio
!
! PROCEDURE:
!       To convert mixing ratio in g/kg to parts-per-million, the following
!       is used:
!                                                 MW(Dry Air)
!         ppmv(MOL) = 1000 . Mixing_Ratio(MOL) . -------------
!                                                   MW(MOL)
!
!       where MW(Dry Air) = Average molecular weight of dry air
!             MW(MOL)     = Molecular weight of the gas in question.
!
!       The factor of 1000 derives from the product of the g/kg to g/g
!       scale factor (0.001) and the "parts-per" to "parts-per-million"
!       scale factor (1.0e+06)
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 03-May-2000
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION MR_to_PPMV_scalar( Mixing_Ratio, &  ! Input
                              Molecule_ID,  &  ! Optional Input
                              Message_Log ) &  ! Error messaging
                            RESULT( ppmv )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Mixing_Ratio
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: ppmv
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MR_to_PPMV'
    REAL(fp),     PARAMETER :: SCALE_FACTOR = G_TO_KG * PPV_TO_PPMV
    ! Local variables
    INTEGER :: Id

    ! Setup
    ppmv = -ONE
    IF ( Mixing_Ratio < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input mixing ratio < 0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( PRESENT( Molecule_ID ) ) THEN
      IF ( Molecule_ID < 1 .OR. Molecule_ID > MAX_N_MOLECULAR_SPECIES ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Unrecognised Molecule_ID.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Id = Molecule_ID
    ELSE
      Id = 1  ! Default value is for water vapor
    END IF

    ! Convert mass to volume mixing ratio
    ppmv = SCALE_FACTOR * Mixing_Ratio * MW_DRYAIR / MOLECULAR_WEIGHT( Id )

  END FUNCTION MR_to_PPMV_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION MR_to_PPMV_rank1( Mixing_Ratio, &  ! Input
                             Molecule_ID,  &  ! Input
                             Message_Log ) &  ! Error messaging
                           RESULT( ppmv )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Mixing_Ratio(:)
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: ppmv(SIZE(Mixing_Ratio))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MR_to_PPMV'
    ! Local variables
    INTEGER :: k

    ! Setup
    ppmv = -ONE

    ! Convert mass to volume mixing ratio
    DO k = 1, SIZE(Mixing_Ratio)
      ppmv( k ) = MR_to_PPMV_scalar( Mixing_Ratio( k ), &
                                     Molecule_ID=Molecule_ID, &
                                     Message_Log=Message_Log )
      IF ( ppmv( k ) < ZERO ) RETURN
    END DO

  END FUNCTION MR_to_PPMV_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       PPMV_to_MR
!
! PURPOSE:
!       Function to convert gas concentrations from volume mixing ratio in
!       ppmv to mass mixing ratio in g/kg.
!
! CALLING SEQUENCE:
!       Mixing_Ratio = PPMV_to_MR( ppmv,                    &  ! Input
!                                  Molecule_ID=Molecule_ID, &  ! Input
!                                  Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       ppmv:          Volume mixing ratio of gas.
!                      UNITS:      ppmv
!                      TYPE:       REAL(fp)
!                      DIMENSION:  Scalar or Rank-1 (K x 1)
!                      ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Molecule_ID:   HITRAN molecular designation identifying the
!                      molecule for which the concentration units
!                      conversion is required. If not specified, the
!                      default value is that for water vapor.
!                      Valid values are:
!                        1: H2O       9: SO2      17: HI       25: H2O2
!                        2: CO2      10: NO2      18: ClO      26: C2H2
!                        3: O3       11: NH3      19: OCS      27: C2H6
!                        4: N2O      12: HNO3     20: H2CO     28: PH3
!                        5: CO       13: OH       21: HOCl     29: COF2
!                        6: CH4      14: HF       22: N2       30: SF6
!                        7: O2       15: HCl      23: HCN      31: H2S
!                        8: NO       16: HBr      24: CH3Cl    32: HCOOH
!                      UNITS:      N/A
!                      TYPE:       INTEGER
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: OPTIONAL, INTENT(IN)
!
!       Message_Log:   Character string specifying a filename in which any
!                      Messages will be logged. If not specified, or if an
!                      error occurs opening the log file, the default action
!                      is to output Messages to standard output.
!                      UNITS:      N/A
!                      TYPE:       CHARACTER(*)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! FUNCTION RESULT:
!       Mixing_Ratio:  Mass mixing ratio of gas. If an error occurs,
!                      -1.0 is returned.
!                      UNITS:      g/kg
!                      TYPE:       REAL(fp)
!                      DIMENSION:  Same as input ppmv
!
! PROCEDURE:
!       To convert ppmv to mixing ratio, the following is used:
!
!                                          MW(MOL)
!         mr(MOL) = 0.001 . ppmv(MOL) . -------------
!                                        MW(Dry Air)
!
!       where MW(Dry Air) = Average molecular weight of dry air
!             MW(MOL)     = Molecular weight of the gas in question.
!
!       The factor of 0.001 derives from the product of the g/g to g/kg
!       scale factor (1000) and the "parts-per-million" to "parts-per"
!       scale factor (1.0e-06)
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 03-May-2000
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION PPMV_to_MR_scalar( ppmv,         &  ! Input
                              Molecule_ID,  &  ! Optional Input
                              Message_Log ) &  ! Optional input
                            RESULT( Mixing_Ratio )
    ! Arguments
    REAL(fp),               INTENT(IN) :: ppmv
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Mixing_Ratio
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PPMV_to_MR'
    REAL(fp),     PARAMETER :: SCALE_FACTOR = KG_TO_G * PPMV_TO_PPV
    ! Local variables
    INTEGER :: Id

    ! Setup
    Mixing_Ratio = -ONE
    IF ( ppmv < ZERO ) THEN
      Mixing_Ratio = -ONE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input ppmv < 0.0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( PRESENT( Molecule_ID ) ) THEN
      IF ( Molecule_ID < 1 .OR. Molecule_ID > MAX_N_MOLECULAR_SPECIES ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Unrecognised Molecule_ID.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Id = Molecule_ID
    ELSE
      Id = 1  ! Default value is for water vapor
    END IF

    ! Convert volume to mass mixing ratio
    Mixing_Ratio = SCALE_FACTOR * ppmv * MOLECULAR_WEIGHT( Id ) / MW_DRYAIR

  END FUNCTION PPMV_to_MR_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION PPMV_to_MR_rank1( ppmv,         &  ! Input
                             Molecule_ID,  &  ! Optional Input
                             Message_Log ) &  ! Error messaging
                           RESULT( Mixing_Ratio )
    ! Arguments
    REAL(fp),               INTENT(IN) :: ppmv(:)
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Mixing_Ratio(SIZE(ppmv))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PPMV_to_MR'
    ! Local variables
    INTEGER :: k

    ! Setup
    Mixing_Ratio = -ONE

    ! Convert volume to mass mixing ratio
    DO k = 1, SIZE(ppmv)
      Mixing_Ratio( k ) = PPMV_to_MR_scalar( ppmv( k ), &
                                             Molecule_ID=Molecule_ID, &
                                             Message_Log=Message_Log )
      IF ( Mixing_Ratio( k ) < ZERO ) RETURN
    END DO

  END FUNCTION PPMV_to_MR_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       PPMV_to_PP
!
! PURPOSE:
!       Function to convert gas concentrations from ppmv to partial pressure
!
! CALLING SEQUENCE:
!       Partial_Pressure = PPMV_to_PP( Pressure,                &  ! Input
!                                      ppmv,                    &  ! Input
!                                      Water_Vapor=Water_Vapor, &  ! Optional Input
!                                      Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Total atmospheric pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!
!       ppmv:              Gas volume mixing ratio in ppmv.
!                          UNITS:      ppmv
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!   
! OPTIONAL INPUT ARGUMENTS:
!       Water_Vapor:       Water vapor partial pressure. If this argument is
!                          not supplied, the mandatory input PPMV argument is
!                          assumed to be water vapor.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!                          ATTRIBUTES: INTENT(IN)
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
!       Partial_Pressure:  Gas partial pressure. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       To convert volume mixing ratio in ppmv of a molecular species, 
!       designated by MOL, to partial pressure, the following is used,
!
!         pp(MOL) = 1.0e-06 . ppmv(MOL) . ( Pressure - pp(H2O) )
!
!       If the input molecule is water vapor, the partial pressure is
!       determined using,
!                                                                 1
!         pp(H2O) = 1.0e-06 . ppmv(H2O) . Pressure . -----------------------------
!                                                     1 + ( 1.0e-06 . ppmv(H20) )
!
!       This is based on the same methods used in the LBLRTM conversion code
!       in its lblatm.f module where molecular mixing ratios (either ppmv or
!       g/kg) are always with respect to DRY AIR.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 03-May-2000
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION PPMV_to_PP_scalar( Pressure,     &  ! Input
                              ppmv,         &  ! Input
                              Water_Vapor,  &  ! Optional Input
                              Message_Log ) &  ! Error messaging
                            RESULT( Partial_Pressure )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: ppmv
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Partial_Pressure
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PPMV_to_PP'
    ! Local variables
    REAL(fp) :: ppv
    REAL(fp) :: Dry_Air_Pressure

    ! Setup
    Partial_Pressure = -ONE
    IF ( Pressure < TOLERANCE .OR. &
         ppmv     < ZERO           ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input Pressure = 0.0 or ppmv < or = 0.0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Convert input to parts-per
    ppv = PPMV_TO_PPV * ppmv

    ! Calculate the dry air partial pressure
    IF ( PRESENT(Water_Vapor) ) THEN
      IF ( Water_Vapor < ZERO ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Input Water_Vapor partial pressure < 0.0.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Dry_Air_Pressure = Pressure - Water_Vapor
    ELSE
      Dry_Air_Pressure = Pressure * ( ONE / ( ONE + ppv ) )
    END IF

    ! Calculate the partial pressure
    Partial_Pressure = ppv * Dry_Air_Pressure

  END FUNCTION PPMV_to_PP_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION PPMV_to_PP_rank1( Pressure,     &  ! Input
                             ppmv,         &  ! Input
                             Water_Vapor,  &  ! Optional Input
                             Message_Log ) &  ! Error messaging
                           RESULT( Partial_Pressure )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: ppmv(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Partial_Pressure(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PPMV_to_PP'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    Partial_Pressure = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(ppmv) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Pressure/ppmv array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
      
    ! Calculate the partial pressure
    IF ( PRESENT(Water_Vapor) ) THEN
    
      ! Water vapour argument is present
      IF ( SIZE(Water_Vapor) /= n ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Inconsistent input Pressure/Water_Vapor array sizes.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      ENDIF
      DO k = 1, n
        Partial_Pressure( k ) = PPMV_to_PP_scalar( Pressure( k ), &
                                                   ppmv( k ), &
                                                   Water_Vapor=Water_Vapor( k ), &
                                                   Message_Log=Message_Log )
        IF ( Partial_Pressure( k ) < ZERO ) RETURN
      END DO

    ELSE
    
      ! Water vapour argument is NOT present
      DO k = 1, n
        Partial_Pressure( k ) = PPMV_to_PP_scalar( Pressure( k ), &
                                                   ppmv( k ), &
                                                   Message_Log=Message_Log )
        IF ( Partial_Pressure( k ) < ZERO ) RETURN
      END DO

    END IF

  END FUNCTION PPMV_to_PP_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       PP_to_PPMV
!
! PURPOSE:
!       Function to convert gas concentrations from partial pressure to 
!       volume mixing ratio in ppmv
!
! CALLING SEQUENCE:
!       ppmv = PP_to_PPMV( Pressure,                &  ! Input
!                          Partial_Pressure,        &  ! Input
!                          Water_Vapor=Water_Vapor, &  ! Optional Input
!                          Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Total atmospheric pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!
!       Partial_Pressure:  Gas partial pressure
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as pressure
!                          ATTRIBUTES: INTENT(IN)
!   
! OPTIONAL INPUT ARGUMENTS:
!       Water_Vapor:       Water vapor volume mixing ratio. If this argument is
!                          not supplied, the mandatory Partial_Pressure argument
!                          is assumed to be water vapor.
!                          UNITS:      ppmv
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
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
!       ppmv:              Gas volume mixing ratio. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      ppmv
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       To convert the partial pressure of a molecular species, designated
!       by MOL, to volume mixing ratio in ppmv, the following is used,
!
!                                     pp(MOL)
!         ppmv(MOL) = 1.0e+06 . ------------------
!                                Dry_Air_Pressure
!
!       where
!                                                     1
!         Dry_Air_Pressure = Pressure . ----------------------------
!                                        1 + ( ppmv(H2O) . 1.0e-6 )
!
!       If the input molecule is water vapor, the dry air pressure is
!       determined simply using,
!
!         Dry_Air_Pressure = Pressure - pp(H2O)
!
!       This is based on the same methods used in the LBLRTM conversion code
!       in its lblatm.f module where molecular mixing ratios (either ppmv or
!       g/kg) are always with respect to DRY AIR.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 03-May-2000
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION PP_to_PPMV_scalar( Pressure,         &  ! Input
                              Partial_Pressure, &  ! Input
                              Water_Vapor,      &  ! Optional Input
                              Message_Log )     &  ! Error messaging
                            RESULT( ppmv )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Partial_Pressure
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: ppmv
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PP_to_PPMV'
    ! Local variables
    REAL(fp) :: ppv
    REAL(fp) :: Dry_Air_Pressure

    ! Setup
    ppmv = -ONE
    IF ( Pressure         < TOLERANCE .OR. &
         Partial_Pressure < ZERO           ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input pressures < or = 0.0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate the dry air partial pressure
    IF ( PRESENT(Water_Vapor) ) THEN
      IF ( Water_Vapor < ZERO ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Input Water_Vapor ppmv < 0.0.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      ! Convert water vapor ppmv to ppv
      ppv = PPMV_TO_PPV * Water_Vapor
      Dry_Air_Pressure = Pressure * ( ONE / ( ONE + ppv ) )
    ELSE
      ! Input partial pressure is for water vapor
      Dry_Air_Pressure = Pressure - Partial_Pressure
    END IF

    ! Calculate the volume mixing ratio
    ppmv = PPV_TO_PPMV * Partial_Pressure / Dry_Air_Pressure

  END FUNCTION PP_to_PPMV_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION PP_to_PPMV_rank1( Pressure,         &  ! Input
                             Partial_Pressure, &  ! Input
                             Water_Vapor,      &  ! Optional Input
                             Message_Log )     &  ! Error messaging
                           RESULT( ppmv )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Partial_Pressure(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: ppmv(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PP_to_PPMV'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    ppmv = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Partial_Pressure) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Pressure/Partial_Pressure array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate the volume mixing ratio
    IF ( PRESENT(Water_Vapor) ) THEN

      ! Water vapour argument is present
      IF ( SIZE(Water_Vapor) /= n ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Inconsistent input Pressure/Water_Vapor array sizes.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      ENDIF
      DO k = 1, n
        ppmv( k ) = PP_to_PPMV_scalar( Pressure( k ), &
                                       Partial_Pressure( k ), &
                                       Water_Vapor=Water_Vapor( k ), &
                                       Message_Log=Message_Log )
        IF ( ppmv( k ) < ZERO ) RETURN
      END DO

    ELSE

      ! Water vapour argument is NOT present
      DO k = 1, n
        ppmv( k ) = PP_to_PPMV_scalar( Pressure( k ), &
                                       Partial_Pressure( k ), &
                                       Message_Log=Message_Log )
        IF ( ppmv( k ) < ZERO ) RETURN
      END DO

    END IF

  END FUNCTION PP_to_PPMV_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       MR_to_PP
!
! PURPOSE:
!       Function to convert gas concentrations from mixing ratio in g/kg
!       to partial pressure in hectoPascals.
!
! CALLING SEQUENCE:
!       Partial_Pressure = MR_to_PP( Pressure,                &  ! Input
!                                    Mixing_Ratio,            &  ! Input
!                                    Molecule_ID=Molecule_ID, &  ! Optional Input
!                                    Water_Vapor=Water_Vapor, &  ! Optional Input
!                                    Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Total atmospheric pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!
!       Mixing_Ratio:      Mass mixing ratio of gas.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as pressure
!                          ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Molecule_ID:       HITRAN molecular designation identifying the
!                          molecule for which the concentration units
!                          conversion is required. If not specified, the
!                          default value is that for water vapor.
!                          Valid values are:
!                            1: H2O       9: SO2      17: HI       25: H2O2
!                            2: CO2      10: NO2      18: ClO      26: C2H2
!                            3: O3       11: NH3      19: OCS      27: C2H6
!                            4: N2O      12: HNO3     20: H2CO     28: PH3
!                            5: CO       13: OH       21: HOCl     29: COF2
!                            6: CH4      14: HF       22: N2       30: SF6
!                            7: O2       15: HCl      23: HCN      31: H2S
!                            8: NO       16: HBr      24: CH3Cl    32: HCOOH
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!
!       Water_Vapor:       Water vapor partial pressure. If this argument is
!                          not supplied, the mandatory MIXING_RATIO argument
!                          is assumed to be water vapor.
!                          This argument is ignored if the specified or default
!                          molecule ID is set to 1.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
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
!       Partial_Pressure:  Gas partial pressure. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       First a definition. In this routine, the mass mixing ratio of a gas
!       is defined to be, for a given volume, the ratio of the mass of the gas 
!       in question to that of DRY air.
!
!       The ideal gas equation can be written as,
! 
!         p = n.R0.T     ................................................(1)
!
!       for unit volume, where 
!         n  = number of moles of gas,
!         R0 = universal gas constant
!         T  = temperature.
!
!       The mass mixing ratio of a gas is defined as the mass of the gas with
!       respect to the mass of dry air in the same volume. If we use eqn(1)
!       to construct expressions for the partial pressures of a particular
!       gas and dry air, we get,
!
!         pp(MOL)    = n(MOL).R0.T     .................................(2a)
!
!       and
!
!         pp(DryAir) = n(DryAir).R0.T     ..............................(2b)
!
!
!       Dividing eqn(2a) by (2b) and rearranging we get,
!
!                     n(MOL) 
!         pp(MOL) = ----------- . pp(DryAir)     ........................(3)
!                    n(DryAir)
!
!       Replacing the expssion for the number of moles of a substance into
!       eqn(3) gives us,
!
!                     m(MOL)     MW(DryAir)
!         pp(MOL) = --------- . ------------ . pp(DryAir)
!                    MW(MOL)      m(DryAir)
!
!                     m(MOL)       MW(DryAir)
!                 = ----------- . ------------ . pp(DryAir)
!                    m(DryAir)      MW(MOL)
!
!                                     MW(DryAir)
!                 = 0.001 . w(MOL) . ------------ . pp(DryAir)     ......(4)
!                                      MW(MOL)
!
!       where m(MOL)     = mass of gas MOL in grams,
!             m(DryAir)  = mass of dry air in grams,
!             MW(MOL)    = molecular weight of of gas MOL in grams,
!             MW(DryAir) = effective molecular weight of dry air in grams,
!             w(MOL)     = mass mixing ratio of gas MOL in g/kg.
!
!       The factor of 0.001 in eqn(4) is to convert the units of the mixing
!       ratio from g/kg to g/g.
!
!       Two cases need to be addressed:
!
!
!       -------------
!       1) MOL == H2O
!       -------------
!
!       If the gas for which the mixing ratio is to be converted is water
!       vapor, then eqn(4) can be written as,
!
!                   
!         pp(H2O) = WX . ( p(Total) - pp(H2O) )     .....................(5)
!
!                                    MW(DryAir)
!       where WX = 0.001 . w(H2O) . ------------
!                                     MW(H2O)
!
!       Expanding eqn(5) further, we get
!                   
!         pp(H2O) = ( WX . p(Total) ) - ( WX . pp(H2O) )
!
!       and,
!
!         pp(H2O) . ( 1 + WX ) = WX . p(Total)
!
!       and finally,
!
!                        WX
!         pp(H2O) = ------------ . p(Total)     .........................(6)
!                    ( 1 + WX )
!
!       Eqn(6) is used to determine the partial pressure for water vapor
!       input in this routine.
!
!
!       -------------------
!       2) MOL is *not* H2O
!       -------------------
!
!       Using eqn(4) to determine an expression for the non-water vapor
!       gases, we get,
!
!         pp(MOL) = WX . p(DryAir)     ..................................(7)
!
!                                    MW(DryAir)
!       where WX = 0.001 . w(MOL) . ------------
!                                     MW(MOL)
!
!       and w(MOL) is still defined as in eqn(4).
!
!
!       This is based on the same methods used in the LBLRTM conversion code
!       in its lblatm.f module where molecular mixing ratios (either ppmv or
!       g/kg) are always with respect to DRY AIR.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 20-Sep-2002
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION MR_to_PP_scalar( Pressure,     &  ! Input
                            Mixing_Ratio, &  ! Input
                            Molecule_ID,  &  ! Optional Input
                            Water_Vapor,  &  ! Optional Input
                            Message_Log ) &  ! Error messaging
                          RESULT( Partial_Pressure )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Mixing_Ratio
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Partial_Pressure
    ! Local parameters
    CHARACTER(*),  PARAMETER :: ROUTINE_NAME = 'MR_to_PP'
    ! Local variables
    INTEGER :: Id
    REAL(fp) :: w

    ! Setup
    Partial_Pressure = -ONE
    IF ( Pressure     < TOLERANCE .OR. &
         Mixing_Ratio < ZERO      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input Pressure/Mixing_Ratio < or = 0.0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( PRESENT( Molecule_ID ) ) THEN
      IF ( Molecule_ID < 1 .OR. Molecule_ID > MAX_N_MOLECULAR_SPECIES ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Unrecognised Molecule_ID.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Id = Molecule_ID
    ELSE
      Id = 1  ! Default value is for water vapor
    END IF

    ! Calculate the "w" factor
    w = G_TO_KG * Mixing_Ratio * MW_DRYAIR / MOLECULAR_WEIGHT( Id )

    ! Convert amount based on molecule ID
    IF ( Id == 1 ) THEN

      ! Conversion is for water vapor
      Partial_Pressure = ( w / ( ONE + w ) ) * Pressure

    ELSE

      ! Conversion is for some other molecular species
      IF ( .NOT. PRESENT(Water_Vapor) ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Must specify Water_Vapor partial pressure for '//&
                              'non-H2O Mixing_Ratio input.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      IF ( Water_Vapor < ZERO ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Input Water_Vapor partial pressure < 0.0.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Partial_Pressure = w * ( Pressure - Water_Vapor )

    END IF

  END FUNCTION MR_to_PP_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION MR_to_PP_rank1( Pressure,     &  ! Input
                           Mixing_Ratio, &  ! Input
                           Molecule_ID,  &  ! Optional Input
                           Water_Vapor,  &  ! Opotional Input
                           Message_Log ) &  ! Error messaging
                         RESULT( Partial_Pressure )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Mixing_Ratio(:)
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Partial_Pressure(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MR_to_PP'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    Partial_Pressure = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Mixing_Ratio) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Pressure/Mixing_Ratio array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate partial pressures
    IF ( PRESENT(Water_Vapor) ) THEN

      ! Water vapour argument is present
      IF ( SIZE(Water_Vapor) /= n ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Inconsistent input Pressure/Water_Vapor array sizes.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      ENDIF
      DO k = 1, n
        Partial_Pressure( k ) = MR_to_PP_scalar( Pressure( k ), &
                                                 Mixing_Ratio( k ), &
                                                 Molecule_ID=Molecule_ID, &
                                                 Water_Vapor=Water_Vapor( k ), &
                                                 Message_Log=Message_Log )
        IF ( Partial_Pressure( k ) < ZERO ) RETURN
      END DO

    ELSE

      ! Water vapor argument was NOT passed
      DO k = 1, n
        Partial_Pressure( k ) = MR_to_PP_scalar( Pressure( k ), &
                                                 Mixing_Ratio( k ), &
                                                 Molecule_ID=Molecule_ID, &
                                                 Message_Log=Message_Log )
        IF ( Partial_Pressure( k ) < ZERO ) RETURN
      END DO

    END IF

  END FUNCTION MR_to_PP_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       PP_to_MR
!
! PURPOSE:
!       Function to convert gas concentrations from partial pressure in
!       hectoPascals to mass mixing ratio in g/kg.
!
! CALLING SEQUENCE:
!       Mixing_Ratio = PP_to_MR( Pressure,                &  ! Input
!                                Partial_Pressure,        &  ! Input
!                                Molecule_ID=Molecule_ID, &  ! Optional Input
!                                Water_Vapor=Water_Vapor, &  ! Optional Input
!                                Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Total atmospheric pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!
!       Partial_Pressure:  Partial pressure of gas.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Molecule_ID:       HITRAN molecular designation identifying the
!                          molecule for which the concentration units
!                          conversion is required. If not specified, the
!                          default value is that for water vapor.
!                          Valid values are:
!                            1: H2O       9: SO2      17: HI       25: H2O2
!                            2: CO2      10: NO2      18: ClO      26: C2H2
!                            3: O3       11: NH3      19: OCS      27: C2H6
!                            4: N2O      12: HNO3     20: H2CO     28: PH3
!                            5: CO       13: OH       21: HOCl     29: COF2
!                            6: CH4      14: HF       22: N2       30: SF6
!                            7: O2       15: HCl      23: HCN      31: H2S
!                            8: NO       16: HBr      24: CH3Cl    32: HCOOH
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!
!       Water_Vapor:       Water vapor mass mixing ratio. If this argument is
!                          not supplied, the mandatory PARTIAL_PRESSURE argument
!                          is assumed to be water vapor.
!                          This argument is ignored if the specified or default
!                          molecule ID is set to 1.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!
!       Message_Log:       Character string specifying a filename in which any
!                          messages will be logged. If not specified, or if an
!                          error occurs opening the log file, the default action
!                          is to output messages to standard output.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! FUNCTION RESULT:
!       Mixing_Ratio:      Mass mixing ratio of gas. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       First a definition. In this routine, the mass mixing ratio of a gas
!       is defined to be, for a given volume, the ratio of the mass of the gas 
!       in question to that of DRY air.
!
!       The ideal gas equation can be written as,
! 
!         p = n.R0.T     ................................................(1)
!
!       for unit volume, where 
!         n  = number of moles of gas,
!         R0 = universal gas constant
!         T  = temperature.
!
!       The mass mixing ratio of a gas is defined as the mass of the gas with
!       respect to the mass of dry air in the same volume. If we use eqn(1)
!       to construct expressions for the partial pressures of a particular
!       gas and dry air, we get,
!
!         pp(MOL)    = n(MOL).R0.T     .................................(2a)
!
!       and
!
!         pp(DryAir) = n(DryAir).R0.T     ..............................(2b)
!
!
!       Dividing eqn(2a) by (2b) and rearranging we get,
!
!                     n(MOL) 
!         pp(MOL) = ----------- . pp(DryAir)     ........................(3)
!                    n(DryAir)
!
!       Replacing the expssion for the number of moles of a substance into
!       eqn(3) gives us,
!
!                     m(MOL)     MW(DryAir)
!         pp(MOL) = --------- . ------------ . pp(DryAir)
!                    MW(MOL)      m(DryAir)
!
!                     m(MOL)       MW(DryAir)
!                 = ----------- . ------------ . pp(DryAir)
!                    m(DryAir)      MW(MOL)
!
!                                     MW(DryAir)
!                 = 0.001 . w(MOL) . ------------ . pp(DryAir)     ......(4)
!                                      MW(MOL)
!
!       where m(MOL)     = mass of gas MOL in grams,
!             m(DryAir)  = mass of dry air in grams,
!             MW(MOL)    = molecular weight of of gas MOL in grams,
!             MW(DryAir) = effective molecular weight of dry air in grams,
!             w(MOL)     = mass mixing ratio of gas MOL in g/kg.
!
!       The factor of 0.001 in eqn(4) is to convert the units of the mixing
!       ratio from g/kg to g/g.
!
!       Thus to determine the mixing ratio of the gas in question from its
!       partial pressure, eqn(4) is rearranged to give,
!
!                            MW(MOL)        pp(MOL)
!         w(MOL) = 1000 . ------------ . ------------     ...............(5)
!                          MW(DryAir)     pp(DryAir)
!
!
!       Two cases need to be addressed:
!
!
!       -------------
!       1) MOL == H2O
!       -------------
!
!       If the gas for which the partial pressure is to be converted is water
!       vapor, then the dry air partial pressure required in eqn(5) is simply
!       computed using,
!
!         pp(DryAir) = p(Total) - pp(H2O)
!
!       which is then used in eqn(5) to compute the water vapour mixing ratio.
!
!
!       -------------------
!       2) MOL is *not* H2O
!       -------------------
!
!       For this gas, the dry air partial pressure must be computed using the 
!       water vapor mixing ratio supplied in the optional Water_Vapor argument.
!
!       Eqn(4) above can be rewritten as,
!
!         pp(H2O) = WX . p(DryAir)     ..................................(6)
!
!                                    MW(DryAir)
!       where WX = 0.001 . w(H2O) . ------------
!                                     MW(H2O)
!
!       Rearranging eqn(6), we then get,
!
!                      pp(H2O)
!         p(DryAir) = ---------     .....................................(7)
!                        WX
!
!       But, eqn(6) can also be written as,
!
!         pp(H2O) = WX . ( p(Total) - pp(H2O) )
!
!                 = ( WX . p(Total) ) - ( WX . pp(H2O) )
!
!       and thus,
!
!         pp(H2O) . ( 1 + WX ) = WX . p(Total)
!
!       with finally,
!
!                        WX
!         pp(H2O) = ------------ . p(Total)     .........................(8)
!                    ( 1 + WX )
!
!       Substituting eqn(8) into eqn(7) gives,
!
!                      p(Total)
!         p(DryAir) = ----------     ....................................(9)
!                      1 + WX 
!
!       Eqn(9) is used to compute the dry air partial pressure which is then
!       used in enq(5) to compute the mixing ratio of the gas in question.
!
!
!       This is based on the same methods used in the LBLRTM conversion code
!       in its lblatm.f module where molecular mixing ratios (either ppmv or
!       g/kg) are always with respect to DRY AIR.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 20-Sep-2002
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION PP_to_MR_scalar( Pressure,         &  ! Input
                            Partial_Pressure, &  ! Input
                            Molecule_ID,      &  ! Optional Input
                            Water_Vapor,      &  ! Optional Input
                            Message_Log )     &  ! Error messaging
                          RESULT( Mixing_Ratio )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Partial_Pressure
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Mixing_Ratio
    ! Local parameters
    CHARACTER(*),  PARAMETER :: ROUTINE_NAME = 'PP_to_MR'
    ! Local variables
    INTEGER :: Id
    REAL(fp) :: w
    REAL(fp) :: Dry_Air_Pressure

    ! Setup
    Mixing_Ratio = -ONE
    IF ( Pressure         < TOLERANCE .OR. &
         Partial_Pressure < ZERO      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input Pressure/Partial_Pressure < or = 0.0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( PRESENT( Molecule_ID ) ) THEN
      IF ( Molecule_ID < 1 .OR. Molecule_ID > MAX_N_MOLECULAR_SPECIES ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Unrecognised Molecule_ID.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Id = Molecule_ID
    ELSE
      Id = 1  ! Default value is for water vapor
    END IF



    ! Calculate the dry air partial pressure
    IF ( Id == 1 ) THEN

      ! Conversion is for water vapor
      Dry_Air_Pressure = Pressure - Partial_Pressure

    ELSE

      ! Conversion is for some other molecular species
      IF ( .NOT. PRESENT( Water_Vapor ) ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Must specify Water_Vapor mixing ratio for '//&
                              'non-H2O Partial_Pressure input.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      IF ( Water_Vapor < ZERO ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Input Water_Vapor mixing ratio < 0.0.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      w = G_TO_KG * Water_Vapor * MW_DRYAIR / MW_H2O
      Dry_Air_Pressure = Pressure / ( ONE + w ) 

    END IF

    ! Calculate the mass mixing ratio
    w = KG_TO_G * MOLECULAR_WEIGHT( Id ) / MW_DRYAIR
    Mixing_Ratio = w * Partial_Pressure / Dry_Air_Pressure

  END FUNCTION PP_to_MR_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION PP_to_MR_rank1( Pressure,     &  ! Input
                           Mixing_Ratio, &  ! Input
                           Molecule_ID,  &  ! Optional Input
                           Water_Vapor,  &  ! Opotional Input
                           Message_Log ) &  ! Error messaging
                         RESULT( Partial_Pressure )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Mixing_Ratio(:)
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Partial_Pressure(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PP_to_MR'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    Partial_Pressure = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Mixing_Ratio) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Pressure/Mixing_Ratio array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate partial pressures
    IF ( PRESENT(Water_Vapor) ) THEN

      ! Water vapour argument is present
      IF ( SIZE(Water_Vapor) /= n ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Inconsistent input Pressure/Water_Vapor array sizes.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      ENDIF
      DO k = 1, n
        Partial_Pressure( k ) = PP_to_MR_scalar( Pressure( k ), &
                                                 Mixing_Ratio( k ), &
                                                 Molecule_ID=Molecule_ID, &
                                                 Water_Vapor=Water_Vapor( k ), &
                                                 Message_Log=Message_Log )
        IF ( Partial_Pressure( k ) < ZERO ) RETURN
      END DO

    ELSE

      ! Water vapor argument was NOT passed
      DO k = 1, n
        Partial_Pressure( k ) = PP_to_MR_scalar( Pressure( k ), &
                                                 Mixing_Ratio( k ), &
                                                 Molecule_ID=Molecule_ID, &
                                                 Message_Log=Message_Log )
        IF ( Partial_Pressure( k ) < ZERO ) RETURN
      END DO

    END IF

  END FUNCTION PP_to_MR_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       PP_to_MD
!
! PURPOSE:
!       Function to convert gas concentrations in pressure units to 
!       mass density.
!
! CALLING SEQUENCE:
!       Mass_Density = PP_to_MD( Pressure,                &  ! Input
!                                Temperature,             &  ! Input
!                                Molecule_ID=Molecule_ID, &  ! Optional Input
!                                Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Total or partial pressure to provide number
!                          density of air or specific gas species.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1
!                          ATTRIBUTES: INTENT(IN)
!   
!       Temperature:       Atmospheric temperature
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as pressure
!                          ATTRIBUTES: INTENT(IN)
!   
! OPTIONAL INPUT ARGUMENTS:
!       Molecule_ID:       HITRAN molecular designation identifying the
!                          molecule for which the concentration units
!                          conversion is required. If not specified, the
!                          default value is that for water vapor.
!                          Valid values are:
!                            1: H2O       9: SO2      17: HI       25: H2O2
!                            2: CO2      10: NO2      18: ClO      26: C2H2
!                            3: O3       11: NH3      19: OCS      27: C2H6
!                            4: N2O      12: HNO3     20: H2CO     28: PH3
!                            5: CO       13: OH       21: HOCl     29: COF2
!                            6: CH4      14: HF       22: N2       30: SF6
!                            7: O2       15: HCl      23: HCN      31: H2S
!                            8: NO       16: HBr      24: CH3Cl    32: HCOOH
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
!       Mass_Density:      Number density of gas. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      g/m^3
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       The ideal gas law is
!
!         p.V = n.R.T
!
!       where,
!         p = Pressure
!         V = volume
!         n = number of moles of the gas,
!         R = molar gas constant,
!         T = Temperature.
!
!       This can be written in terms of the number of molecules, N,
!
!                 N
!         p.V = ---- . R.T                  ..............(1)
!                NA
!
!       where NA = Avogadro's constant (number of molecules in 1 mole.)
!
!       Eqn.(1) can be recast to provide the number of molecules in the 
!       volume, the number density, nd,
!
!               N     p.NA
!         nd = --- = ------  molecules/m^3  ..............(2)
!               V     R.T
!
!       The mass density is related to the number density by the simple
!       relation,
!
!                    MW
!         md = nd . ----                    ..............(3)
!                    NA
!
!       Substituting equation (2) into (3) gives,
!
!               p.MW
!         md = ------  g/m^3                ..............(4)
!               R.T
!
!       The value of p0 used in this routine is expressed in pascals (Pa) whereas
!       the input partial Pressure is expected to be hPa (more common unit). Thus
!       there is a factor of 100 to include,
!
!               100.p.MW
!         md = ----------  g/m^3
!                 R.T
!
! UNITS ANALYSIS:
!       p  :  kg.m^-1.s^-2  (Pascals)
!       MW :  g.mol^-1
!       R  :  J.mol^-1.K^-1 == kg.m^2.s^-2.mol^-1.K^-1
!       T  :  K
!
!                  kg       g      s^2.mol.K     1
!         md == ------- . ----- . ----------- . ---
!                m.s^2     mol      kg.m^2       K
!
!                 g
!            == -----
!                m^3
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 08-Nov-2004
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION PP_to_MD_scalar( Pressure,     &  ! Input
                            Temperature,  &  ! Input
                            Molecule_ID,  &  ! Optional Input
                            Message_Log ) &  ! Error messaging
                          RESULT( Mass_Density )
    ! Arguments
    REAL(fp),          INTENT(IN) :: Pressure
    REAL(fp),          INTENT(IN) :: Temperature
    INTEGER,        OPTIONAL, INTENT(IN) :: Molecule_ID
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Mass_Density
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PP_to_MD'
    ! Local variables
    INTEGER :: Id

    ! Setup
    Mass_Density = -ONE
    IF ( Pressure    < ZERO      .OR. &
         Temperature < TOLERANCE      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input partial pressure < 0, or temperature = 0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( PRESENT( Molecule_ID ) ) THEN
      IF ( Molecule_ID < 1 .OR. Molecule_ID > MAX_N_MOLECULAR_SPECIES ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Unrecognised Molecule_ID.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Id = Molecule_ID
    ELSE
      Id = 1  ! Default value is for water vapor
    END IF

    ! Convert partial pressure to mass density
    Mass_Density = HPA_TO_PA * Pressure * MOLECULAR_WEIGHT( Id ) / ( Temperature * R0 )

  END FUNCTION PP_to_MD_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION PP_to_MD_rank1( Pressure,     &  ! Input
                           Temperature,  &  ! Input
                           Molecule_ID,  &  ! Optional Input
                           Message_Log ) &  ! Error messaging
                         RESULT( Mass_Density )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Mass_Density(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PP_to_MD'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    Mass_Density = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Pressure/Temperature array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Convert partial pressure to mass density
    DO k = 1, n
      Mass_Density( k ) = PP_to_MD_scalar( Pressure( k ), &
                                           Temperature( k ), &
                                           Molecule_ID=Molecule_ID, &
                                           Message_Log=Message_Log )
      IF ( Mass_Density( k ) < ZERO ) RETURN
    END DO

  END FUNCTION PP_to_MD_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       MD_to_PP
!
! PURPOSE:
!       Function to convert gas concentration mass density to partial
!       pressure.
!
! CALLING SEQUENCE:
!       Pressure = MD_to_PP( Mass_Density,            &  ! Input
!                            Temperature,             &  ! Input
!                            Molecule_ID=Molecule_ID, &  ! Optional Input
!                            Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Mass_Density:      Number density of gas.
!                          UNITS:      g/m^3
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1
!                          ATTRIBUTES: INTENT(IN)
!   
!       Temperature:       Atmospheric temperature
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Mass_Density
!                          ATTRIBUTES: INTENT(IN)
!   
! OPTIONAL INPUT ARGUMENTS:
!       Molecule_ID:       HITRAN molecular designation identifying the
!                          molecule for which the concentration units
!                          conversion is required. If not specified, the
!                          default value is that for water vapor.
!                          Valid values are:
!                            1: H2O       9: SO2      17: HI       25: H2O2
!                            2: CO2      10: NO2      18: ClO      26: C2H2
!                            3: O3       11: NH3      19: OCS      27: C2H6
!                            4: N2O      12: HNO3     20: H2CO     28: PH3
!                            5: CO       13: OH       21: HOCl     29: COF2
!                            6: CH4      14: HF       22: N2       30: SF6
!                            7: O2       15: HCl      23: HCN      31: H2S
!                            8: NO       16: HBr      24: CH3Cl    32: HCOOH
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
!       Pressure:          Partial pressure for specified gas species.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Mass_Density
!
! PROCEDURE:
!       The ideal gas law is
!
!         p.V = n.R.T
!
!       where,
!         p = Pressure
!         V = volume
!         n = number of moles of the gas,
!         R = molar gas constant,
!         T = Temperature.
!
!       This can be written in terms of the number of molecules, N,
!
!                 N
!         p.V = ---- . R.T                  ..............(1)
!                NA
!
!       where NA = Avogadro's constant (number of molecules in 1 mole.)
!
!       Eqn.(1) can be recast to provide the number of molecules in the 
!       volume, the number density, nd,
!
!               N     p.NA
!         nd = --- = ------  molecules/m^3  ..............(2)
!               V     R.T
!
!       The mass density is related to the number density by the simple
!       relation,
!
!                    MW
!         md = nd . ----                    ..............(3)
!                    NA
!
!       Substituting equation (2) into (3) gives,
!
!               p.MW
!         md = ------  g/m^3                ..............(4)
!               R.T
!
!       and inverting equation (4) gives,
!
!              md.R.T
!         p = -------- Pa
!                MW
!
!       The pressure result above is determined in units of Pascals, so there
!       is a factor of 100 to include to return the pressure in units of
!       hectoPascals,
!
!              md.R.T
!         p = -------- hPa
!              100.MW
!
! UNITS ANALYSIS:
!       md :  g.m^-3
!       R  :  J.mol^-1.K^-1 == kg.m^2.s^-2.mol^-1.K^-1
!       T  :  K
!       MW :  g.mol^-1
!
!                g        kg.m^2          mol
!         p == ----- . ----------- . K . -----
!               m^3     s^2.mol.K          g
!
!                 kg
!           == -------
!               m.s^2
!
!               kg.m.s^-2      F
!           == ----------- == --- == pressure.
!                  m^2         A
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 08-Nov-2004
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION MD_to_PP_scalar( Mass_Density, &  ! Input
                            Temperature,  &  ! Input
                            Molecule_ID,  &  ! Optional Input
                            Message_Log ) &  ! Error messaging
                          RESULT( Pressure )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Mass_Density
    REAL(fp),               INTENT(IN) :: Temperature
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Pressure
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MD_to_PP'
    ! Local variables
    INTEGER :: Id

    ! Setup
    Pressure = -ONE
    IF ( Mass_Density < ZERO      .OR. &
         Temperature  < TOLERANCE      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input mass density < 0, or temperature = 0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( PRESENT( Molecule_ID ) ) THEN
      IF ( Molecule_ID < 1 .OR. Molecule_ID > MAX_N_MOLECULAR_SPECIES ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Unrecognised Molecule_ID.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Id = Molecule_ID
    ELSE
      Id = 1  ! Default value is for water vapor
    END IF

    ! Convert mass density to partial pressure
    Pressure = PA_TO_HPA * Mass_Density * R0 * Temperature / MOLECULAR_WEIGHT( Id )

  END FUNCTION MD_to_PP_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION MD_to_PP_rank1( Mass_Density, &  ! Input
                           Temperature,  &  ! Input
                           Molecule_ID,  &  ! Optional Input
                           Message_Log ) &  ! Error messaging
                         RESULT( Pressure )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Mass_Density(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    INTEGER,      OPTIONAL, INTENT(IN) :: Molecule_ID
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Pressure(SIZE(Mass_Density))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MD_to_PP'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    Pressure = -ONE
    n = SIZE(Mass_Density)
    IF ( SIZE(Temperature) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Mass Density/Temperature array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    
    ! Convert mass density to partial pressure
    DO k = 1, n
      Pressure( k ) = MD_to_PP_scalar( Mass_Density( k ), &
                                       Temperature( k ), &
                                       Molecule_ID=Molecule_ID, &
                                       Message_Log=Message_Log )
      IF ( Pressure( k ) < ZERO ) RETURN
    END DO

  END FUNCTION MD_to_PP_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       PP_to_ND
!
! PURPOSE:
!       Function to convert gas concentrations in pressure units to 
!       number density.
!
! CALLING SEQUENCE:
!       Number_Density = PP_to_ND( Pressure,                 &  ! Input
!                                  Temperature,              &  ! Input
!                                  Message_Log=Message_Log )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Total or partial pressure to provide number
!                          density of air or specific gas species.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!   
!       Temperature:       Atmospheric temperature
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
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
!       Number_Density:    Number density of gas. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      molecules/m^3
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       The ideal gas law is
!
!         p.V = n.R.T
!
!       where,
!         p = Pressure
!         V = volume
!         n = number of moles of the gas,
!         R = molar gas constant,
!         T = Temperature.
!
!       This can be written in terms of the number of molecules, N,
!
!                 N
!         p.V = ---- . R.T                  ..............(1)
!                NA
!
!       where NA = Avogadro's constant (number of molecules in 1 mole.)
!
!       Eqn.(1) can be recast to provide the number of molecules in the 
!       volume, the number density, nd,
!
!               N     p.NA
!         nd = --- = ------  molecules/m^3  ..............(2)
!               V     R.T
!
!       At standard Temperature and Pressure (T0=273.15K, p0=101325Pa), this
!       number density is know as the Loschmidt constant, L0, the molecular
!       density of 1 mole of an ideal gas. Thus we have the generic form of
!       eqn.(2) and the STP form,
!
!               p0.NA
!         L0 = -------  molecules/m^3       ..............(3)
!               R.T0
!
!       Taking the ratio of eqns.(2) and (3) gives,
!       
!         nd    p.NA     R.T0
!         -- = ------ . -------
!         L0    R.T      p0.NA
!
!       and rearranging gives,
!
!                    p      T0
!         nd = L0 . ---- . ----  molecules/m^3
!                    p0     T 
!
!       The value of p0 used in this routine is expressed in pascals (Pa) whereas
!       the input partial Pressure is expected to be hPa (more common unit). Thus
!       there is a factor of 100 to include,
!
!                    100.p     T0
!         nd = L0 . ------- . ----  molecules/m^3
!                     p0       T 
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 15-Nov-2001
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION PP_to_ND_scalar( Pressure,     &  ! Input
                            Temperature,  &  ! Input
                            Message_Log ) &  ! Error messaging
                          RESULT( Number_Density )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Temperature
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Number_Density
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PP_to_ND'

    ! Setup
    Number_Density = -ONE
    IF ( Pressure    < ZERO      .OR. &
         Temperature < TOLERANCE      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input partial Pressure < 0, or Temperature = 0.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Convert partial pressure to number density
    Number_Density = HPA_TO_PA * Pressure * L0 * T0 / ( Temperature * P0 )

  END FUNCTION PP_to_ND_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION PP_to_ND_rank1( Pressure,     &  ! Input
                           Temperature,  &  ! Input
                           Message_Log ) &  ! Error messaging
                         RESULT( Number_Density )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Number_Density(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PP_to_ND'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    Number_Density = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Pressure/Temperature array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Convert partial pressure to number density
    DO k = 1, n
      Number_Density( k ) = PP_to_ND_scalar( Pressure( k ), &
                                             Temperature( k ), &
                                             Message_Log=Message_Log )
      IF ( Number_Density( k ) < ZERO ) RETURN
    END DO

  END FUNCTION PP_to_ND_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       ND_to_PP
!
! PURPOSE:
!       Function to convert gas concentrations from number density to
!       pressure units.
!
! CALLING SEQUENCE:
!       Partial_Pressure = ND_to_PP( Number_Density,         &  ! Input
!                                    Temperature,            &  ! Input
!                                    Message_Log=Message_Log )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Number_Density:    Molecular density.
!                          UNITS:      molecules/m^3
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!   
!       Temperature:       Atmospheric Temperature
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Number_Density
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
!       Partial_Pressure:  Gas partial pressure. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Number_Density
!
! PROCEDURE:
!       The ideal gas law is
!
!         p.V = n.R.T
!
!       where,
!         p = Pressure
!         V = volume
!         n = number of moles of the gas,
!         R = molar gas constant,
!         T = Temperature.
!
!       This can be written in terms of the number of molecules, N,
!
!                 N
!         p.V = ---- . R.T                  ..............(1)
!                NA
!
!       where NA = Avogadro's constant (number of molecules in 1 mole.)
!
!       Recasting eqn.(1) to provide the Pressure,
!
!              nd
!         p = ---- . R.T  Pa                ..............(2)
!              NA
!
!       where nd = the number density in molecules/m^3 = N/V.
!
!       At standard Temperature and Pressure (T0=273.15K, p0=101325Pa), the
!       number density of eqn.(2) is known as the Loschmidt constant, L0,
!       the molecular density of 1 mole of an ideal gas. Thus we have the
!       generic form of eqn.(2) and the STP form,
!
!               L0
!         p0 = ---- . R.T0  Pa              ..............(3)
!               NA
!
!       Taking the ratio of eqns.(2) and (3) gives,
!       
!          p      nd     T
!         ---- = ---- . ----  Pa
!          p0     L0     T0
!
!       and rearranging gives,
!
!                   nd      T
!         p = p0 . ---- . ----  Pa
!                   L0     T0 
!
!       The value of p0 used in this routine is expressed in pascals (Pa) whereas
!       the output Pressure is returned as hPa (more common unit). Thus there
!       is a factor of 100 to include,
!
!                          nd      T
!         p = 0.01 . p0 . ---- . ----  hPa
!                          L0     T0 
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 15-Nov-2001
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION ND_to_PP_scalar( Number_Density, &  ! Input
                            Temperature,    &  ! Input
                            Message_Log )   &  ! Optional input
                          RESULT( Pressure )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Number_Density
    REAL(fp),               INTENT(IN) :: Temperature
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Pressure
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ND_to_PP'

    ! Setup
    Pressure = -ONE
    IF ( Number_Density < ZERO      .OR. &
         Temperature    < TOLERANCE      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input number density < 0.0, or Temperature < or = 0.0', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Convert number density to pressure
    Pressure = PA_TO_HPA * P0 * Number_Density * Temperature / ( L0 * T0 )

  END FUNCTION ND_to_PP_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION ND_to_PP_rank1( Number_Density, &  ! Input
                           Temperature,    &  ! Input
                           Message_Log )   &  ! Optional input
                         RESULT( Pressure )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Number_Density(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Pressure(SIZE(Number_Density))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ND_to_PP'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    Pressure = -ONE
    n = SIZE(Number_Density)
    IF ( SIZE(Temperature) /= n ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Number_Density/Temperature array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Convert number density to pressure
    DO k = 1, n
      Pressure( k ) = ND_to_PP_scalar( Number_Density( k ), &
                                       Temperature( k ), &
                                       Message_Log=Message_Log )
      IF ( Pressure( k ) < ZERO ) RETURN
    END DO

  END FUNCTION ND_to_PP_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       PPMV_to_ND
!
! PURPOSE:
!       Function to convert volume mixing ratio in ppmv to number density.
!
! CALLING SEQUENCE:
!       Number_Density = PPMV_to_ND( Pressure,                &  ! Input
!                                    Temperature,             &  ! Input
!                                    ppmv,                    &  ! Input
!                                    Water_Vapor=Water_Vapor, &  ! Optional Input
!                                    Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:         Total atmospheric pressure.
!                         UNITS:      hectoPascals, hPa
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar or Rank-1 (K x 1)
!                         ATTRIBUTES: INTENT(IN)
!
!       Temperature:      Atmospheric temperature
!                         UNITS:      Kelvin, K
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Same as Pressure
!                         ATTRIBUTES: INTENT(IN)
!
!       ppmv:             Volume mixing ratio.
!                         UNITS:      ppmv
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Same as Pressure
!                         ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Water_Vapor:      Water vapor number density. If this argument is
!                         not supplied, the mandatory PPMV input argument
!                         is assumed to be water vapor.
!                         UNITS:      molecules/m^3
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar or Rank-1 (K x 1)
!                         ATTRIBUTES: INTENT(IN)
!
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
!       Number_Density:   Number density of gas. If an error occurs,
!                         -1.0 is returned.
!                         UNITS:      molecules/m^3
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       To convert the volume mixing ratio of a molecular species, designated
!       by MOL, to number density, the following is used,
!
!         nd(MOL) = 1.0e-06 . ppmv(MOL) . ( nd(TOT) - nd(H2O) )     .....(1)
!
!       where
!                         Pressure          T0
!         nd(TOT) = L0 . ---------- . -------------  molecules/m^3
!                            p0        Temperature 
!
!       and L0 = Loschmidt number,
!           p0 = Standard pressure,
!           T0 = Standard temperature.
!
!       If the input molecule is water vapor, the number density is
!       determined using,
!
!                        ppmv(H2O) . 1.0e-06
!         nd(H2O) = ----------------------------- . nd(TOT)     .....(2)
!                    1 + ( ppmv(H2O) . 1.0e-06 )
!
!       Rearranging eqn.(2) gives the same form as eqn.(1).
!
!       This is based on the same methods used in the LBLRTM conversion code
!       in its lblatm.f module where molecular mixing ratios (either ppmv or
!       g/kg) are always with respect to DRY AIR.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 15-Nov-2001
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION PPMV_to_ND_scalar( Pressure,     &  ! Input
                              Temperature,  &  ! Input
                              ppmv,         &  ! Input
                              Water_Vapor,  &  ! Optional input
                              Message_Log ) &  ! Error messaging
                            RESULT( Number_Density )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Temperature
    REAL(fp),               INTENT(IN) :: ppmv
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Number_Density
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PPMV_to_ND'
    ! Local variables
    REAL(fp) :: ppv
    REAL(fp) :: Total_Density

    ! Setup
    Number_Density = -ONE
    IF ( Pressure    < ZERO      .OR. &
         Temperature < TOLERANCE .OR. &
         ppmv        < ZERO           ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input pressure,ppmv < 0.0, or temperature < or = 0.0', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    END IF

    ! Convert ppmv to ppv
    ppv = PPMV_TO_PPV * ppmv

    ! Calculate total air number density
    Total_Density = PP_to_ND_scalar( Pressure, &
                                     Temperature, &
                                     Message_Log=Message_Log )
    IF ( Total_Density < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error calculating total number density.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    END IF

    ! Calculate the molecular number density
    IF ( PRESENT(Water_Vapor) ) THEN

      ! Water vapor argument is present
      IF ( Water_Vapor < ZERO ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Input water vapor number density < 0.0', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      Number_Density = ppv * ( Total_Density - Water_Vapor )

    ELSE

      ! Water vapor argument was NOT passed
      Number_Density = ( ppv / ( ONE + ppv ) ) * Total_Density

    END IF

  END FUNCTION PPMV_to_ND_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION PPMV_to_ND_rank1( Pressure,     &  ! Input
                             Temperature,  &  ! Input
                             ppmv,         &  ! Input
                             Water_Vapor,  &  ! Optional input
                             Message_Log ) &  ! Error messaging
                           RESULT( Number_Density )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),               INTENT(IN) :: ppmv(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: Number_Density(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PPMV_to_ND'
    ! Local variables
    INTEGER :: k, n

    ! Setup
    Number_Density = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature) /= n .OR. &
         SIZE(ppmv)        /= n      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate the molecular number density
    IF ( PRESENT(Water_Vapor) ) THEN

      ! Water vapor argument was passed
      IF ( SIZE(Water_Vapor) /= n ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Inconsistent input Pressure/Water_Vapor array sizes.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      ENDIF
      DO k = 1, n
        Number_Density( k ) = PPMV_to_ND_scalar( Pressure( k ), &
                                                 Temperature( k ), &
                                                 ppmv( k ), &
                                                 Water_Vapor=Water_Vapor( k ), &
                                                 Message_Log=Message_Log )
        IF ( Number_Density( k ) < ZERO ) RETURN
      END DO

    ELSE

      ! Water vapor argument was *NOT* passed
      DO k = 1, n
        Number_Density( k ) = PPMV_to_ND_scalar( Pressure( k ), &
                                                 Temperature( k ), &
                                                 ppmv( k ), &
                                                 Message_Log=Message_Log )
        IF ( Number_Density( k ) < ZERO ) RETURN
      END DO

    END IF

  END FUNCTION PPMV_to_ND_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       ND_to_PPMV
!
! PURPOSE:
!       Function to convert gas concentrations from number density to volume
!       mixing ratio in ppmv.
!
! CALLING SEQUENCE:
!       ppmv = ND_to_PPMV( Pressure,                &  ! Input
!                          Temperature,             &  ! Input
!                          Number_Density,          &  ! Input
!                          Water_Vapor=Water_Vapor, &  ! Optional Input
!                          Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Total atmospheric pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!   
!       Temperature:       Atmospheric temperature
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!   
!       Number_Density:    Molecular number density.
!                          UNITS:      molecules/m^3
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!   
! OPTIONAL INPUT ARGUMENTS:
!       Water_Vapor:       Water vapor volume mixing ratio. If this argument is
!                          not supplied, the mandatory NUMBER_DENSITY input
!                          argument is assumed to be water vapor.
!                          UNITS:      ppmv
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
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
!       ppmv:              Gas volume mixing ratio. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      ppmv
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       To convert the number density of a molecular species, designated
!       by MOL, to volume mixing ratio, the following is used,
!
!                                  nd(MOL)
!         ppmv(MOL) = 1.0e+06 . -------------
!                                nd(DRY_AIR)
!
!       where
!
!         nd(DRY_AIR) = nd(TOT) - nd(H2O)
!
!       and
!                         Pressure          T0
!         nd(TOT) = L0 . ---------- . -------------  molecules/m^3
!                            p0        Temperature 
!
!       with L0 = Loschmidt number,
!            p0 = Standard pressure,
!            T0 = Standard temperature, and
!
!       If the input molecule is NOT water vapor, the dry air number density
!       is determined using,
!
!                                     1
!         nd(DRY_AIR) = ----------------------------- . nd(TOT)
!                        1 + ( ppmv(H2O) . 1.0e-06 )
!
!       This is based on the same methods used in the LBLRTM conversion code
!       in its lblatm.f module where molecular mixing ratios (either ppmv or
!       g/kg) are always with respect to DRY AIR.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 15-Nov-2001
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION ND_to_PPMV_scalar( Pressure,       &  ! Input
                              Temperature,    &  ! Input
                              Number_Density, &  ! Input
                              Water_Vapor,    &  ! Optional Input
                              Message_Log )   &  ! Error messaging
                            RESULT( ppmv )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Temperature
    REAL(fp),               INTENT(IN) :: Number_Density
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: ppmv
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ND_to_PPMV'
    ! Local variables
    REAL(fp) :: Total_Density
    REAL(fp) :: ppv
    REAL(fp) :: Dry_Air_Density
    
    ! Setup
    ppmv = -ONE
    IF ( Pressure       < ZERO      .OR. &
         Temperature    < TOLERANCE .OR. &
         Number_Density < ZERO           ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input pressure, number density < 0.0, or Temperature < or = 0.0', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    END IF

    ! Calculate total air number density
    Total_Density = PP_to_ND_scalar( Pressure, &
                                     Temperature, &
                                     Message_Log=Message_Log )
    IF ( Total_Density < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error calculating total number density.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    END IF

    ! Calculate the dry air number density
    IF ( PRESENT(Water_Vapor) ) THEN

      ! Water vapor argument is present
      IF ( Water_Vapor < ZERO ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Input water vapor volume mixing ratio < 0.0.', &
                              FAILURE, &
                             Message_Log=Message_Log )
        RETURN
      END IF
      ppv = PPMV_TO_PPV * Water_Vapor
      Dry_Air_Density = Total_Density * ( ONE / ( ONE + ppv ) )

    ELSE

      ! Water vapor argument was *NOT* passed
      Dry_Air_Density = Total_Density - Number_Density

    END IF

    ! Calculate the volume mixing ratio
    ppmv = PPV_TO_PPMV * Number_Density / Dry_Air_Density

  END FUNCTION ND_to_PPMV_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION ND_to_PPMV_rank1( Pressure,       &  ! Input
                             Temperature,    &  ! Input
                             Number_Density, &  ! Input
                             Water_Vapor,    &  ! Optional Input
                             Message_Log )   &  ! Error messaging
                           RESULT( ppmv )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),               INTENT(IN) :: Number_Density(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: ppmv(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'ND_to_PPMV'
    ! Local variables
    INTEGER :: k, n
    
    ! Setup
    ppmv = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature)    /= n .OR. &
         SIZE(Number_Density) /= n      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input Pressure/Number_Density/Temperature array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate the volume mixing ratio
    IF ( PRESENT(Water_Vapor) ) THEN

      ! Water vapor argument was passed
      IF ( SIZE(Water_Vapor) /= n ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Inconsistent input Pressure/Water_Vapor array sizes.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      DO k = 1, n
        ppmv( k ) = ND_to_PPMV_scalar( Pressure( k ), &
                                       Temperature( k ), &
                                       Number_Density( k ), &
                                       Water_Vapor=Water_Vapor( k ), &
                                       Message_Log=Message_Log )
        IF ( ppmv( k ) < ZERO ) RETURN
      END DO

    ELSE

      ! Water vapor argument was *NOT* passed
      DO k = 1, n
        ppmv( k ) = ND_to_PPMV_scalar( Pressure( k ), &
                                       Temperature( k ), &
                                       Number_Density( k ), &
                                       Message_Log=Message_Log )
        IF ( ppmv( k ) < ZERO ) RETURN
      END DO

    END IF

  END FUNCTION ND_to_PPMV_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       PPMV_to_KMOL
!
! PURPOSE:
!       Function to convert layer gas concentrations from volume mixing ratio
!       in ppmv to column density in kmol.cm^-2.
!
! CALLING SEQUENCE:
!       kmol_per_cm2 = PPMV_to_KMOL( Pressure,                &  ! Input
!                                    Temperature,             &  ! Input
!                                    Delta_Z,                 &  ! Input
!                                    ppmv,                    &  ! Input
!                                    Water_Vapor=Water_Vapor, &  ! Optional Input
!                                    Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Average layer pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!
!       Temperature:       Average layer temperature
!                          UNITS:      Kelvin, K (K)
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!
!       Delta_Z:           Layer thickness
!                          UNITS:      metres, m
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!
!       ppmv:              Average layer gas concentration
!                          UNITS:      ppmv
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!   
! OPTIONAL INPUT ARGUMENTS:
!       Water_Vapor:       Water vapor column density. If this argument is
!                          not supplied, the mandatory PPMV input
!                          argument is assumed to be water vapor.
!                          UNITS:      kmol.cm^-2
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
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
!       kmol_per_cm2:      Layer gas concentration. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      kmol.cm^-2
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same input Pressure
!
! PROCEDURE:
!       The column density of a particular molecular species, designated
!       by MOL, given its volume mixing ratio is then given by,
!
!         cd(MOL) = 1.0e-06 . ppmv(MOL) . cd(DRY_AIR)  kmol.cm^-2
!
!       First, the total number density is calculated using,
!
!                         p      T0
!         nd(TOT) = L0 . ---- . ----  molecules.m^-3
!                         p0     T 
!
!       and L0 = Loschmidt number,
!           p0 = Standard pressure,
!           T0 = Standard temperature.
!
!       The total column density is calculated by multiplying the number
!       density by the layer thickness in metres, dz,
!
!         cd(TOT) = dz . nd(TOT)  molecules.m^-2
!
!       This result is scaled by:
!         a) 10^-4 to convert molecules.m^-2 to molecules.cm^-2
!         b) 1/Na to convert molecules.cm^-2 to mol.cm^-2, and
!         c) 10^-3 to convert mol.cm^-2 to kmol.cm^-2
!       giving,
!                    1.0e-07
!         cd(TOT) = --------- . dz . nd(TOT)  kmol.cm^-2
!                      Na
!      
!       The dry air column density is then calculated using,
!
!                                    1
!         cd(DRY_AIR) = ----------------------------- . cd(TOT)  kmol.cm^-2
!                        1 + ( ppmv(H2O) . 1.0e-06 )
!
!       for water vapor ppmv input or,
!
!         cd(DRY_AIR) = cd(TOT) - cd(H2O)
!
!       if cd(H2O) is supplied via the optional Water_Vapor argument.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 04-Jan-2001
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION PPMV_to_KMOL_scalar( Pressure,     &  ! Input
                                Temperature,  &  ! Input
                                Delta_Z,      &  ! Input
                                ppmv,         &  ! Input
                                Water_Vapor,  &  ! Optional input
                                Message_Log ) &  ! Error messaging
                              RESULT( kmol_per_cm2 )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Temperature
    REAL(fp),               INTENT(IN) :: Delta_Z
    REAL(fp),               INTENT(IN) :: ppmv
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: kmol_per_cm2
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PPMV_to_KMOL'
    REAL(fp),     PARAMETER :: SCALE_FACTOR = 1.0e-07_fp / NA
    ! Local variables
    REAL(fp) :: ppv
    REAL(fp) :: Number_Density
    REAL(fp) :: Column_Density
    
    ! Setup
    kmol_per_cm2 = -ONE
    IF ( Pressure     < TOLERANCE .OR. &
         Temperature  < TOLERANCE .OR. &
         ABS(Delta_Z) < TOLERANCE .OR. &
         ppmv         < ZERO           ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input values < or = 0.0 found.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Convert ppmv to ppv
    ppv = PPMV_TO_PPV * ppmv

    ! Calculate total air number density in molecules.m^-3
    Number_Density = PP_to_ND_scalar( Pressure, &
                                      Temperature, &
                                      Message_Log=Message_Log )
    IF ( Number_Density < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error calculating total number density.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    END IF

    ! Calculate total air column density in kmol.cm^-2
    Column_Density = SCALE_FACTOR * ABS(Delta_Z) * Number_Density

    ! Calculate the gas column density
    IF ( PRESENT( Water_Vapor ) ) THEN
      IF ( Water_Vapor < ZERO ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Input water vapor column density < 0.0', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      kmol_per_cm2 = ppv * ( Column_Density - Water_Vapor )
    ELSE
      kmol_per_cm2 = ( ppv / ( ONE + ppv ) ) * Column_Density
    END IF

  END FUNCTION PPMV_to_KMOL_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION PPMV_to_KMOL_rank1( Pressure,     &  ! Input
                               Temperature,  &  ! Input
                               Delta_Z,      &  ! Input
                               ppmv,         &  ! Input
                               Water_Vapor,  &  ! Optional input
                               Message_Log ) &  ! Error messaging
                             RESULT( kmol_per_cm2 )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),               INTENT(IN) :: Delta_Z(:)
    REAL(fp),               INTENT(IN) :: ppmv(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: kmol_per_cm2(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'PPMV_to_KMOL'
    ! Local variables
    INTEGER :: k, n
    
    ! Setup
    kmol_per_cm2 = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature) /= n .OR. & 
         SIZE(Delta_Z)     /= n .OR. & 
         SIZE(ppmv)        /= n      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    END IF

    ! Calculate the gas column density
    IF ( PRESENT(Water_Vapor) ) THEN

      ! Water vapor argument was passed
      IF ( SIZE(Water_Vapor) /= n ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Inconsistent input Pressure/Water_Vapor array sizes.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      DO k = 1, n
        kmol_per_cm2( k ) = PPMV_to_KMOL_scalar( Pressure( k ), &
                                                 Temperature( k ), &
                                                 Delta_Z( k ), &
                                                 ppmv( k ), &
                                                 Water_Vapor=Water_Vapor( k ), &
                                                 Message_Log=Message_Log )
        IF ( kmol_per_cm2( k ) < ZERO ) RETURN
      END DO

    ELSE

      ! Water vapor argument was *NOT*passed
      DO k = 1, n
        kmol_per_cm2( k ) = PPMV_to_KMOL_scalar( Pressure( k ), &
                                                 Temperature( k ), &
                                                 Delta_Z( k ), &
                                                 ppmv( k ), &
                                                 Message_Log=Message_Log )
        IF ( kmol_per_cm2( k ) < ZERO ) RETURN
      END DO

    END IF
   
  END FUNCTION PPMV_to_KMOL_rank1


!------------------------------------------------------------------------------
!S+
! NAME:
!       KMOL_to_PPMV
!
! PURPOSE:
!       Function to convert layer gas concentrations from column density in
!       kmol.cm^-2 to volume mixing ratio in ppmv.
!
! CALLING SEQUENCE:
!       ppmv = KMOL_to_PPMV( Pressure,                &  ! Input
!                            Temperature,             &  ! Input
!                            Delta_Z,                 &  ! Input
!                            kmol_per_cm2,            &  ! Input
!                            Water_Vapor=Water_Vapor, &  ! Optional Input
!                            Message_Log=Message_Log  )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Pressure:          Average layer pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or Rank-1 (K x 1)
!                          ATTRIBUTES: INTENT(IN)
!
!       Temperature:       Average layer temperature
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!
!       Delta_Z:           Layer thickness
!                          UNITS:      metres, m
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!
!       kmol_per_cm2:      Column density for molecular species.
!                          UNITS:      kmol.cm^-2
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!   
! OPTIONAL INPUT ARGUMENTS:
!       Water_Vapor:       Water vapor volume mixing ratio. If this argument
!                          is not supplied, the mandatory KMOL_PER_CM2 input
!                          argument is assumed to be water vapor.
!                          UNITS:      ppmv
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as pressure
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
!       ppmv:              Layer gas volume mixing ratio. If an error occurs,
!                          -1.0 is returned.
!                          UNITS:      ppmv
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as input Pressure
!
! PROCEDURE:
!       The volume mixing ratio of a particular molecular species, designated
!       by MOL, given it's column density is given by,
!
!                                  cd(MOL)
!         ppmv(MOL) = 1.0e-06 . -------------
!                                cd(DRY_AIR)
!
!       First, the total number density is calculated using,
!
!                         p      T0
!         nd(TOT) = L0 . ---- . ----  molecules.m^-3
!                         p0     T 
!
!       and L0 = Loschmidt number,
!           p0 = Standard pressure,
!           T0 = Standard temperature.
!
!       The total column density is calculated by multiplying the number
!       density by the layer thickness in metres, dz,
!
!         cd(TOT) = dz . nd(TOT)  molecules.m^-2
!
!       This result is scaled by:
!         a) 10^-4 to convert molecules.m^-2 to molecules.cm^-2
!         b) 1/Na to convert molecules.cm^-2 to mol.cm^-2, and
!         c) 10^-3 to convert mol.cm^-2 to kmol.cm^-2
!       giving,
!                    1.0e-07
!         cd(TOT) = --------- . dz . nd(TOT)  kmol.cm^-2
!                      Na
!      
!       The dry air column density is then calculated using,
!
!                                    1
!         cd(DRY_AIR) = ----------------------------- . cd(TOT)  kmol.cm^-2
!                        1 + ( ppmv(H2O) . 1.0e-06 )
!
!       for water vapor ppmv input or,
!
!         cd(DRY_AIR) = cd(TOT) - cd(H2O)
!
!       if cd(H2O) is supplied via the optional Water_Vapor argument.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 04-Jan-2001
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  ! ==============
  ! Scalar version
  ! ==============
  FUNCTION KMOL_to_PPMV_scalar( Pressure,     &  ! Input
                                Temperature,  &  ! Input
                                Delta_Z,      &  ! Input
                                kmol_per_cm2, &  ! Input
                                Water_Vapor,  &  ! Optional Input
                                Message_Log ) &  ! Error messaging
                              RESULT( ppmv )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure
    REAL(fp),               INTENT(IN) :: Temperature
    REAL(fp),               INTENT(IN) :: Delta_Z
    REAL(fp),               INTENT(IN) :: kmol_per_cm2
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: ppmv
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'KMOL_to_PPMV'
    REAL(fp),     PARAMETER :: SCALE_FACTOR = 1.0e-07_fp / NA
    ! Local variables
    REAL(fp) :: ppv
    REAL(fp) :: Number_Density
    REAL(fp) :: Column_density
    REAL(fp) :: Dry_Column_density
    
    ! Setup
    ppmv = -ONE
    IF ( Pressure     < TOLERANCE .OR. &
         Temperature  < TOLERANCE .OR. &
         ABS(Delta_Z) < TOLERANCE .OR. &
         kmol_per_cm2 < ZERO      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Input values  < or = 0.0 found.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate total air number density in molecules.m^-3
    Number_Density = PP_to_ND_scalar( Pressure, &
                                      Temperature, &
                                      Message_Log=Message_Log )
    IF ( Number_Density < ZERO ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error calculating total number density.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    END IF

    ! Calculate total air column density in kmol.cm^-2
    Column_Density = SCALE_FACTOR * ABS(Delta_Z) * Number_Density

    ! Calculate the dry gas column density
    IF ( PRESENT(Water_Vapor) ) THEN

      ! Water vapor input was specified
      IF ( Water_Vapor < ZERO ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Input water vapor ppmv < 0.0', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      ppv = PPMV_TO_PPV * Water_Vapor
      Dry_Column_Density = ( ONE / ( ONE + ppv ) ) * Column_Density

    ELSE

      ! Water vapor input was *NOT* specified
      Dry_Column_Density = Column_Density - kmol_per_cm2

    END IF

    ! Convert the column density ratio to ppmv
    ppmv = PPV_TO_PPMV * kmol_per_cm2 / Dry_Column_Density

  END FUNCTION KMOL_to_PPMV_scalar


  ! ==============
  ! Rank-1 version
  ! ==============
  FUNCTION KMOL_to_PPMV_rank1( Pressure,     &  ! Input
                               Temperature,  &  ! Input
                               Delta_Z,      &  ! Input
                               kmol_per_cm2, &  ! Input
                               Water_Vapor,  &  ! Optional Input
                               Message_Log ) &  ! Error messaging
                             RESULT( ppmv )
    ! Arguments
    REAL(fp),               INTENT(IN) :: Pressure(:)
    REAL(fp),               INTENT(IN) :: Temperature(:)
    REAL(fp),               INTENT(IN) :: Delta_Z(:)
    REAL(fp),               INTENT(IN) :: kmol_per_cm2(:)
    REAL(fp),     OPTIONAL, INTENT(IN) :: Water_Vapor(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    ! Function result
    REAL(fp) :: ppmv(SIZE(Pressure))
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'KMOL_to_PPMV'
    ! Local variables
    INTEGER :: k, n
    
    ! Setup
    ppmv = -ONE
    n = SIZE(Pressure)
    IF ( SIZE(Temperature)  /= n .OR. & 
         SIZE(Delta_Z)      /= n .OR. & 
         SIZE(kmol_per_cm2) /= n      ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input array sizes.', &
                            FAILURE, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Convert column desnity to ppmv
    IF ( PRESENT(Water_Vapor) ) THEN

      ! The water vapor argument was passed
      IF ( SIZE(Water_Vapor) /= n ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Inconsistent input Pressure/Water_Vapor array sizes.', &
                              FAILURE, &
                              Message_Log=Message_Log )
        RETURN
      END IF
      DO k = 1, n
        ppmv( k ) = KMOL_to_PPMV_scalar( Pressure( k ), &
                                         Temperature( k ), &
                                         Delta_Z( k ), &
                                         kmol_per_cm2( k ), &
                                         Water_Vapor=Water_Vapor( k ), &
                                         Message_Log=Message_Log )
        IF ( ppmv( k ) < ZERO ) RETURN
      END DO

    ELSE

      ! The water vapor argument was *NOT* passed
      DO k = 1, n
        ppmv( k ) = KMOL_to_PPMV_scalar( Pressure( k ), &
                                         Temperature( k ), &
                                         Delta_Z( k ), &
                                         kmol_per_cm2( k ), &
                                         Message_Log=Message_Log )
        IF ( ppmv( k ) < ZERO ) RETURN
      END DO

    END IF

  END FUNCTION KMOL_to_PPMV_rank1

END MODULE Units_Conversion
