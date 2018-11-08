!
! Level_Layer_Conversion
!
! Module containing routines to convert LEVEL atmospheric profile
! quantities to LAYER quantities.
!
! CONTAINS:
!       Effective_Layer_TP:          Function to calculate the effective (or density
!                                    weighted) temperature and pressure for an
!                                    atmospheric layer.
!
!       Create_Sublevels:            Function to create the sublevels required to
!                                    accurately integrate gas amounts within a layer.
!
!       Integrate_Sublevels:         Function to integrate the temperature and
!                                    absorber amounts to produce average layer
!                                    values.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 01-May-2000
!                       paul.vandelst@ssec.wisc.edu
!

MODULE Level_Layer_Conversion

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Modules used
  USE Type_Kinds, ONLY: fp
  USE Message_Handler
  USE Atmospheric_Properties
  USE Units_Conversion
  USE Profile_Utility_Parameters
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: Effective_Layer_TP
  PUBLIC :: Create_Sublevels
  PUBLIC :: Integrate_Sublevels


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PRIVATE, PARAMETER :: MODULE_RCS_ID = &
  '$Id: Level_Layer_Conversion.f90 2018 2008-04-01 20:12:33Z paul.vandelst@noaa.gov $'


CONTAINS


!------------------------------------------------------------------------------
!S+
! NAME:
!       Effective_Layer_TP
!
! PURPOSE:
!       Function to calculate the effective atmospheric layer temperature and
!       pressure by weighting level values with the integrated layer density.
!
! CALLING SEQUENCE:
!       Error_Status = Effective_Layer_TP( Height,                 &  ! Input
!                                          Pressure,               &  ! Input
!                                          Temperature,            &  ! Input
!                                          Water_Vapor_Pressure,   &  ! Input
!                                          Effective_Pressure,     &  ! Output
!                                          Effective_Temperature,  &  ! Output
!                                          RCS_Id     =RCS_Id,     &  ! Revision control
!                                          Message_Log=Message_Log )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Height:                 Heights of the atmospheric levels.
!                               UNITS:      metres, m
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Rank-1 (N x 1, where N>1)
!                               ATTRIBUTES: INTENT(IN)
!
!       Pressure:               Pressure of the atmospheric levels.
!                               UNITS:      hectoPascals, hPa
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Same as Height
!                               ATTRIBUTES: INTENT(IN)
!
!       Temperature:            Temperature of the atmospheric levels.
!                               UNITS:      Kelvin, K
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Same as Height
!                               ATTRIBUTES: INTENT(IN)
!
!       Water_Vapor_Pressure:   Water vapor partial pressure at the atmospheric levels
!                               UNITS:      hectoPascals, hPa
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Same as Height
!                               ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:            Character string specifying a filename in which any
!                               Messages will be logged. If not specified, or if an
!                               error occurs opening the log file, the default action
!                               is to output Messages to standard output.
!                               UNITS:      N/A
!                               TYPE:       CHARACTER(*)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       Effective_Pressure:     Effective layer pressure.
!                               UNITS:      hectoPascals, hPa
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Rank-1 (N-1 x 1)
!                               ATTRIBUTES: INTENT(OUT)
!
!       Effective_Temperature:  Effective layer temperature.
!                               UNITS:      Kelvin, K
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Same as Effective_Pressure
!                               ATTRIBUTES: INTENT(OUT)
!
!
! OPTIONAL OUTPUT ARGUMENTS:
!       RCS_Id:                 Character string containing the Revision Control
!                               System Id field for the module.
!                               UNITS:      N/A
!                               TYPE:       CHARACTER(*)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(OUT), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:           The return value is an integer defining the error status.
!                               The error codes are defined in the ERROR_HANDLER module.
!                               If == SUCCESS the calculation was successful
!                                  == FAILURE an unrecoverable error occurred
!                               UNITS:      N/A
!                               TYPE:       INTEGER
!                               DIMENSION:  Scalar
!
! PROCEDURE:
!       Based on:
!
!       Gallery, W.O., F.X. Kneizys, and S.A. Clough, "Air mass computer
!         program for atmospheric transmittance/radiance calculation: FSCATM",
!         AFGL-TR-83-0065, 9 March 1983.
!
!       The effective pressure and temperature is defined as,
!               __
!              \           
!               >  p.rho.dz
!              /__         
!         _
!         p = -----------------     ..............................................(1)
!                 __
!                \         
!                 >  rho.dz
!                /__       
!
!       and
!
!               __
!              \           
!               >  T.rho.dz
!              /__         
!         _
!         T = -----------------     ..............................................(2)
!                 __
!                \         
!                 >  rho.dz
!                /__
!
!       where dz == layer thickness.
!
!       Note that the denominators of eqns(1) and (2) can also be referred to as the
!       column density.
!
!       The pressure and total density are both assumed to follow an exponential
!       profile with scale heights H_p and H_rho respectively. For a single layer
!       the numerator of eqn(1) can be written as,
!
!          __ k
!         \                H_p.H_rho
!          >  p.rho.dz = ------------- ( p(k-1).rho(k-1) - p(k).rho(k) )     .....(3)
!         /__             H_p + H_rho
!            k-1
!
!       Similarly for the numerator of eqn(2) using the ideal gas law,
!       p = R_air.rho.T, we get
!
!          __ k
!         \                H_p
!          >  T.rho.dz = -------( p(l-1) - p(l) )     ............................(4)
!         /__             R_air
!            k-1
!
!       and the denominator is given by,
!
!          __ k
!         \           
!          >  rho.dz = H_rho ( rho(l-1) - rho(l) )     ...........................(5)
!         /__         
!            k-1
!
!       where the scale heights are defined as,
!
!                -( z(l) - z(l-1 ) )        
!         H_p = ---------------------     ........................................(6)
!                ln( p(l) / p(l-1) )        
!
!       and
!
!                    -( z(l) - z(l-1 ) )
!         H_rho = -------------------------     ..................................(7)
!                  ln( rho(l) / rho(l-1) )
!
!
!       Note that in eqn.(4) the gas constant is that for *air*, not *dry air*. To
!       determine this the effective molecular weight of air (as a function of pressure)
!       must be determined.
!
!       Breaking down the units of the components, 
!
!         units(p)  = hPa
!                   = 100 Pa
!                   = 100 N.m^-2
!                   = 100 kg.m.s^-2.m^-2
!                   = 100 kg.m^-1.s^-2
!
!                          m2
!         units(eqn(3)) = ----( 100 kg.m^-1.s^-2  .  kg.m^-3 )
!                          m 
!
!                       = 100 kg^2.m^-3.s^-2
!
!                          m  . 100 kg.m^-1.s^-2
!         units(eqn(4)) = -----------------------
!                              J.g^-1.K^-1
!
!                          m  . 100 kg.m^-1.s^-2
!                       = ------------------------
!                           kg.m^2.s^-2.g^-1.K^-1
!
!                       = 100 K.g.m^-2
!                       = 0.1 K.kg.m^-2
!
!         units(eqn(5)) = m  .  kg.m^-3
!                       = kg.m^-2  
!                     
!       So the units of the final equations are:
!
!               _     units(eqn(3))
!         units(p) = ---------------
!                     units(eqn(5))  
!
!                          100 kg^2.m^-3.s^-2
!                  = --------------------
!                         kg.m^-2  
!
!                  = 100 kg.m^-1.s^-2
!                  = 100 kg.m.s^-2.m^-2
!                  = 100 N.m^-2
!                  = 100 Pa
!                  = hPa
!
!               _     units(eqn(4))
!         units(T) = ---------------
!                     units(eqn(5))  
!
!                     0.1 K.kg.m^-2
!                  = ---------------
!                       kg.m^-2  
!
!                  = 0.1 K
!
!       So the final temperatures must be multiplied by 0.1 to get units of K.
!
!       Note for the above units breakdown of eqn(4) that the gas constant
!       for air is computed in units of J.g^-1.K^-1, *NOT* the SI units of
!       J.kg^-1.K^-1. This is done solely to save the conversion operation
!       of g->kg for each loop iteration. Otherwise, for a gas constant in
!       units of J.kg^-1.K^-1, the final scaling factor for the effective
!       temperature would be 100, not 0.1.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC, 03-May-2000
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  FUNCTION Effective_Layer_TP( Height,                &  ! Input
                               Pressure,              &  ! Input
                               Temperature,           &  ! Input
                               Water_Vapor_Pressure,  &  ! Input
                               Effective_Pressure,    &  ! Output
                               Effective_Temperature, &  ! Output
                               RCS_Id,                &  ! Revision control
                               Message_Log )          &  ! Error messaging
                             RESULT( Error_Status )
    ! Arguments
    REAL(fp),               INTENT(IN)  :: Height(:)
    REAL(fp),               INTENT(IN)  :: Pressure(:)
    REAL(fp),               INTENT(IN)  :: Temperature(:)
    REAL(fp),               INTENT(IN)  :: Water_Vapor_Pressure(:)
    REAL(fp),               INTENT(OUT) :: Effective_Pressure(:)
    REAL(fp),               INTENT(OUT) :: Effective_Temperature(:)
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: RCS_Id
    CHARACTER(*), OPTIONAL, INTENT(IN)  :: Message_Log
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Effective_Layer_TP'
    REAL(fp),     PARAMETER :: SCALE_FACTOR = 0.1_fp
    ! Local variables
    CHARACTER(256) :: Message
    INTEGER :: n_Levels, n_Layers
    INTEGER :: k
    REAL(fp) :: MWair
    REAL(fp) :: Rair, Rair_km1, layer_Rair
    REAL(fp) :: RHOair, RHOair_km1
    REAL(fp) :: dz
    REAL(fp) :: H_p, H_rho
    REAL(fp) :: Sum_rho, sum_p_RHO, Sum_T_rho

 
    ! Setup
    ! -----
    Error_Status = SUCCESS
    IF ( PRESENT(RCS_Id) ) RCS_Id = MODULE_RCS_ID

    ! Check input array sizes
    n_Levels = SIZE(Height)
    IF ( SIZE(Pressure)             /= n_Levels .OR. & 
         SIZE(Temperature)          /= n_Levels .OR. & 
         SIZE(Water_Vapor_Pressure) /= n_Levels ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input array sizes.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Check output array sizes
    n_Layers = n_Levels - 1
    IF ( SIZE(Effective_Pressure)    < n_Layers .OR. & 
         SIZE(Effective_Temperature) < n_Layers      ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Output arrays to small to hold result.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Check input array values
    IF ( ANY(Pressure < TOLERANCE) ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input Pressures < or = 0 found.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( ANY(Temperature < TOLERANCE) ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input Temperatures < or = 0 found.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( ANY(Water_Vapor_Pressure < ZERO) ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input water vapor partial pressures < 0 found.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF


    ! Calculate near surface level values
    ! -----------------------------------
    ! Molecular weight of air
    MWair = MW_Air( Pressure( 1 ), &
                    Water_Vapor_Pressure( 1 ), &
                    Message_Log=Message_Log )
    IF ( MWair < ZERO ) THEN
      Error_Status = FAILURE
      WRITE( Message, '( "Error calculating MWair at Level 1. Value = ", es13.6 )' ) &
                      MWair
      CALL Display_Message( ROUTINE_NAME, &
                            TRIM(Message), &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate the gas "constant" in J.g^-1.K^-1
    ! Note that the units are *NOT* SI. The scaling of these
    ! units is addressed in the final scale factor for the
    ! computation of the effective temeprature. (The effective
    ! pressure is not affected as it doesn't require this value)
    Rair_km1 = R0 / MWair

    ! Air density
    RHOair_km1 = Density( Pressure( 1 ), &
                          Temperature( 1 ), &
                          MWair, &
                          Message_Log=Message_Log )

    IF ( RHOair_km1 < ZERO ) THEN
      Error_Status = FAILURE
      WRITE( Message, '( "Error calculating RHOair at Level 1. Value = ", es13.6 )' ) &
                      RHOair_km1
      CALL Display_Message( ROUTINE_NAME, &
                            TRIM(Message), &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF


    ! Loop over layers
    ! ----------------
    Layer_Loop: DO k = 1, n_Layers


      ! Calculate current top of layer values
      ! -------------------------------------
      ! MWair at current Level
      MWair = MW_Air( Pressure( k+1 ), &
                      Water_Vapor_Pressure( k+1 ), &
                      Message_Log=Message_Log )
      IF ( MWair < ZERO ) THEN
        Error_Status = FAILURE
        WRITE( Message, '( "Error calculating MWair at Level ", i4, ". Value = ", es13.6 )' ) &
                        k+1, MWair
        CALL Display_Message( ROUTINE_NAME, &
                              TRIM(Message), &
                              Error_Status, &
                              Message_Log=Message_Log )
        RETURN
      ENDIF

      ! Calculate the gas "constant" in J.g^-1.K^-1
      ! Note that the units are *NOT* SI. The scaling of these
      ! units is addressed in the final scale factor for the
      ! computation of the effective temeprature. (The effective
      ! pressure is not affected as it doesn't require this value)
      Rair = R0 / MWair

      ! Air density at current Level
      RHOair = Density( Pressure( k+1 ), &
                        Temperature( k+1 ), &
                        MWair, &
                        Message_Log=Message_Log )
      IF ( RHOair < ZERO ) THEN
        Error_Status = FAILURE
        WRITE( Message, '( "Error calculating RHOair at Level ", i4, ". Value = ", es13.6 )' ) &
                        k+1, RHOair
        CALL Display_Message( ROUTINE_NAME, &
                              TRIM(Message), &
                              Error_Status, &
                              Message_Log=Message_Log )
        RETURN
      ENDIF


      ! Calculate the layer scale heights
      ! ---------------------------------
      ! Calculate layer thicknesses
      dz = Height( k+1 ) - Height( k )

      ! Pressure scale height
      H_p = dz / LOG(Pressure( k+1 ) / Pressure( k ))

      ! Density scale height
      H_rho = dz / LOG(RHOair / RHOair_km1)


      ! Calculate the effective quantities
      ! ----------------------------------
      ! Calculate the density integral
      Sum_rho = H_rho * ( RHOair - RHOair_km1 )

      ! Effective pressure
      Sum_p_rho = ( ( H_p * H_rho ) / ( H_p + H_rho ) ) * &
                  ( ( Pressure( k+1 ) * RHOair ) - ( Pressure( k ) * RHOair_km1 ) )
      Effective_Pressure( k ) = Sum_p_rho / Sum_rho


      ! Calculate the density weighted layer gas "constant"
      layer_Rair = ( ( Rair_km1*RHOair_km1 ) + ( Rair*RHOair ) ) / ( RHOair_km1+RHOair )

      ! Effective temperature
      Sum_T_rho = ( H_p / layer_Rair ) * ( Pressure( k+1 ) - Pressure( k ) )
      Effective_Temperature( k ) = SCALE_FACTOR * Sum_T_rho / Sum_rho


      ! Save top boundary values for use as bottom boundary
      ! values for next layer
      ! ---------------------------------------------------
      Rair_km1   = Rair
      RHOair_km1 = RHOair

    END DO Layer_Loop

  END FUNCTION Effective_Layer_TP


!------------------------------------------------------------------------------
!S+
! NAME:
!       Create_Sublevels
!
! PURPOSE:
!       Function to create the sublevels used to integrate input profiles
!       to obtain average layer quantities. This routine is called before
!       Integrate_Sublevels.
!
!       Adapted from the UMBC INTLEV.F function supplied with the AIRS RTA.
!
! CALLING SEQUENCE:
!       Error_Status = Create_Sublevels( Level_Pressure,         &  ! Input
!                                        Level_Temperature,      &  ! Input
!                                        Level_Absorber,         &  ! Input
!                                        n_Per_Layer,            &  ! Input
!                                        Sublevel_Pressure,      &  ! Output
!                                        Sublevel_Temperature,   &  ! Output
!                                        Sublevel_Absorber,      &  ! Output
!                                        RCS_Id     =RCS_Id,     &  ! Revision control
!                                        Message_Log=Message_Log )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Level_Pressure:        Pressure of the atmospheric levels.
!                              UNITS:      hectoPascals, hPa
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Rank-1 (K)
!                                            K == number of levels
!                              ATTRIBUTES: INTENT(IN)
!
!       Level_Temperature:     Temperature of the atmospheric levels.
!                              UNITS:      Kelvin, K
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Same as Level_Pressure
!                              ATTRIBUTES: INTENT(IN)
!
!       Level_Absorber:        Absorber concentrations at the atmospheric levels
!                              UNITS:      Doesn't matter - as long as they are
!                                          LEVEL specific.
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Rank-2 (K x J)
!                                            K == number of levels
!                                            J == number of absorbers
!                              ATTRIBUTES: INTENT(IN)
!
!       n_Per_Layer:           Number of sublevels to create in each layer.
!                              Value must be > or = 1.
!                              UNITS:      N/A
!                              TYPE:       INTEGER
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN)
!
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:           Character string specifying a filename
!                              in which any Messages will be logged.
!                              If not specified, or if an error occurs
!                              opening the log file, the default action
!                              is to output Messages to standard output.
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       Sublevel_Pressure:     Pressure of the atmospheric sublevels.
!                              UNITS:      hectoPascals, hPa
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Rank-1 (Ks)
!                                            Ks == number of sublevels
!                                               == ( (K-1) * n_Per_Layer ) + 1
!                              ATTRIBUTES: INTENT(IN)
!
!       Sublevel_Temperature:  Temperature of the atmospheric sublevels.
!                              UNITS:      Kelvin, K
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Same as output Sublevel_Pressure argument
!                              ATTRIBUTES: INTENT(IN)
!
!       Sublevel_Absorber:     Absorber concentrations at the atmospheric Levels
!                              UNITS:      Same as input
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Rank-2 (Ks x J)
!                                            Ks == number of sublevels
!                                            J  == number of absorbers
!                              ATTRIBUTES: INTENT(IN)
!
!
! OPTIONAL OUTPUT ARGUMENTS:
!       RCS_Id:                Character string containing the Revision Control
!                              System Id field for the module.
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(OUT), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:          The return value is an integer defining the error
!                              status. The error codes are defined in the
!                              ERROR_HANDLER module.
!                              If == SUCCESS the calculation was successful
!                                 == FAILURE an unrecoverable error occurred
!                              UNITS:      N/A
!                              TYPE:       INTEGER
!                              DIMENSION:  Scalar
!
! PROCEDURE:
!       The assumption is made that temperature and absorber amount vary
!       linearly with ln(p). To improve the quadrature in integrating Level
!       amounts to a layer value, each input layer, k, is split into N(k)
!       sublayers equally spaced in ln(p),
!
!                                 ln(p[k+1]) - ln(p[k]
!         ln(p[n+1]) - ln(p[n] = ---------------------
!                                          N(k)
!
!       given the pressures, p(1) - p(K) of the input Levels.
!
!       Once the sublevels are defined, the level temperatures and absorber
!       amounts are linearly interpolated at the specific number of sublevels
!       and those interpolates are associated with the sublevel pressures.
!
!       The last corresponding level/sublevel is assigned explicitly.
!
!       Currently, N is independent of k. That is, the same number of sublevels
!       are created for each layer.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 19-Jan-2001
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  FUNCTION Create_Sublevels( Level_Pressure,       &  ! Input
                             Level_Temperature,    &  ! Input
                             Level_Absorber,       &  ! Input
                             n_Per_Layer,          &  ! Input
                             Sublevel_Pressure,    &  ! Output
                             Sublevel_Temperature, &  ! Output
                             Sublevel_Absorber,    &  ! Output
                             RCS_Id,               &  ! Revision control
                             Message_Log )         &  ! Error messaging
                           RESULT( Error_Status )
    ! Arguments
    REAL(fp),               INTENT(IN)  :: Level_Pressure(:)
    REAL(fp),               INTENT(IN)  :: Level_Temperature(:)
    REAL(fp),               INTENT(IN)  :: Level_Absorber(:,:)
    INTEGER,                INTENT(IN)  :: n_Per_Layer
    REAL(fp),               INTENT(OUT) :: Sublevel_Pressure(:)
    REAL(fp),               INTENT(OUT) :: Sublevel_Temperature(:)
    REAL(fp),               INTENT(OUT) :: Sublevel_Absorber(:,:)
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: RCS_Id
    CHARACTER(*), OPTIONAL, INTENT(IN)  :: Message_Log
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Create_Sublevels'
    ! Local variables
    INTEGER :: n_Levels
    INTEGER :: n_Sublevels
    INTEGER :: n_Absorbers
    INTEGER :: i       ! Generic loop/index variable
    INTEGER :: j       ! Absorber index variable
    INTEGER :: k       ! Level index variable
    INTEGER :: n1, n2  ! Sublevel indices within *layer* k
    REAL(fp) :: xn_Per_Layer
    REAL(fp) :: dx
    REAL(fp) :: xn(n_Per_Layer)
    REAL(fp) :: Level_ln_Pressure(SIZE(Level_Pressure))
    REAL(fp) :: Sublevel_ln_Pressure(SIZE(Sublevel_Pressure))


    ! Setup
    ! -----
    Error_Status = SUCCESS
    IF ( PRESENT(RCS_Id) ) RCS_Id = MODULE_RCS_ID

    ! Check size of input arrays
    n_Levels = SIZE(Level_Pressure)
    IF ( SIZE(Level_Temperature)    /= n_Levels .OR. &
         SIZE(Level_Absorber,DIM=1) /= n_Levels      ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input array sizes.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Check sublevel multiplier
    IF ( n_Per_Layer < 1 ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input N_PER_LAYER must be > 0.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Check size of output arrays
    n_Sublevels = ( ( n_Levels - 1 ) * n_Per_Layer ) + 1
    IF ( SIZE(Sublevel_Pressure)       < n_Sublevels .OR. &
         SIZE(Sublevel_Temperature)    < n_Sublevels .OR. &
         SIZE(Sublevel_Absorber,DIM=1) < n_Sublevels      ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Output arrays not large enough to hold result.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    n_Absorbers = SIZE( Level_Absorber, DIM = 2 )
    IF ( SIZE(Sublevel_Absorber,DIM=2) < n_Absorbers ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Output Sublevel_Absorber array does not have '//&
                            'enough absorber dimension elements.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    
    ! Check input Pressure array values
    IF ( ANY(Level_Pressure < TOLERANCE) ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input Pressures < or = 0.0 found.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Calculate the log of the input pressure
    ! ---------------------------------------
    ! Don't really need the WHERE due to the 
    ! input pressure check, but just to be sure
    WHERE ( Level_Pressure > ZERO )
      Level_ln_Pressure = LOG( Level_Pressure )
    ELSEWHERE
      Level_ln_Pressure = ZERO
    END WHERE


    ! Interpolate data
    !
    ! Here we assumes that temperature and absorber amount vary linearly with
    ! the natural logarithm of pressure. Because the interpolation is done at
    ! equal intervals, rather than use a function, the interpolation code is
    ! inline. It's simpler.
    ! -----------------------------------------------------------------------
     
    ! Fill the layer index array
    ! xn = [ 0, 1, 2, ..... n_Per_Layer-1 ]
    xn           = (/ ( REAL(i, fp), i = 0, n_Per_Layer - 1 ) /)
    xn_Per_Layer = REAL(n_Per_Layer, fp)


    ! Loop over layers and linearly interpolate
    ! across layer k (between Levels k and k+1)
    ! to n equally spaced Sublevels.
    !
    !        x(k+1) - x(k)
    !   dx = -------------
    !              n
    ! so that
    !
    !   x = x(k) + ( xn*dx )
    !
    ! where xn and x are vectors.
    !
    ! Note that although the temperature and 
    ! absorber amount are linearly interpolated
    ! between levels, the interpolated values 
    ! are associated with the ln(P) interpolated
    ! values. So, the temperature/absorber
    ! interpolation is effectively exponential.
    Layer_Loop: DO k = 1, n_Levels - 1

      ! Sublevel array indices
      n1 = ( ( k-1 ) * n_Per_Layer ) + 1
      n2 = n1 + n_Per_Layer - 1

      ! Interpolate ln(p)
      dx = ( Level_ln_Pressure( k+1 ) - Level_ln_Pressure( k ) ) / xn_Per_Layer
      Sublevel_ln_Pressure( n1:n2 ) = Level_ln_Pressure( k ) + ( xn * dx )

      ! Interpolate T
      dx = ( Level_Temperature( k+1 ) - Level_Temperature( k ) ) / xn_Per_Layer
      Sublevel_Temperature( n1:n2 ) = Level_Temperature( k ) + ( xn * dx )

      ! Interpolate absorber
      Absorber_Loop: DO j = 1, n_Absorbers
        dx = ( Level_Absorber( k+1, j ) - Level_Absorber( k, j ) ) / xn_Per_Layer
        Sublevel_Absorber( n1:n2, j ) = Level_Absorber( k, j ) + ( xn * dx )
      END DO Absorber_Loop

      ! Convert ln(p) -> p
      Sublevel_Pressure( n1:n2 ) = EXP(Sublevel_ln_Pressure( n1:n2 ))

    END DO Layer_Loop


    ! Assign last Level
    ! -----------------
    Sublevel_Pressure( n_Sublevels )    = Level_Pressure( n_Levels )
    Sublevel_Temperature( n_Sublevels ) = Level_Temperature( n_Levels )
    Sublevel_Absorber( n_Sublevels, : ) = Level_Absorber( n_Levels, : )

  END FUNCTION Create_Sublevels


!------------------------------------------------------------------------------
!S+
! NAME:
!       Integrate_Sublevels
!
! PURPOSE:
!       Function to integrate the sublevel values created by Create_Sublevels
!       to provide average layer temperature and absorber amount. Layer
!       pressure is provided by default, i.e. not from integration.
!
!       Adapted from the UMBC INTEG.F function supplied with the AIRS RTA.
!
! CALLING SEQUENCE:
!       Error_Status = Integrate_Sublevels( Sublevel_Height,        &  ! Input
!                                           Sublevel_Pressure,      &  ! Input
!                                           Sublevel_Temperature,   &  ! Input
!                                           Sublevel_Absorber,      &  ! Input
!                                           n_Per_Layer,            &  ! Input
!                                           H2O_J_Index,            &  ! Input
!                                           Layer_Pressure,         &  ! Output
!                                           Layer_Temperature,      &  ! Output
!                                           Layer_Absorber,         &  ! Output
!                                           RCS_Id     =RCS_Id,     &  ! Revision control
!                                           Message_Log=Message_Log )  ! Error messaging
!
! INPUT ARGUMENTS:
!       Sublevel_Height:       Altitude of the atmospheric sublevels.
!                              UNITS:      metres, m
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Rank-1 (Ks)
!                                            Ks == number of sublevels
!                              ATTRIBUTES: INTENT(IN)
!
!       Sublevel_Pressure:     Pressure of the atmospheric sublevels.
!                              UNITS:      hectoPascals, hPa
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Same as input Sublevel_Height argument
!                              ATTRIBUTES: INTENT(IN)
!
!       Sublevel_Temperature:  Temperature of the atmospheric sublevels.
!                              UNITS:      Kelvin, K
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Same as input Sublevel_Height argument
!                              ATTRIBUTES: INTENT(IN)
!
!       Sublevel_Absorber:     Absorber concentrations at the atmospheric
!                              sublevels
!                              UNITS:      ppmv
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Rank-2 (Ks x J)
!                                            Ks == number of Sublevels
!                                            J  == number of absorbers
!                              ATTRIBUTES: INTENT(IN)
!
!       n_Per_Layer:           Number of sublevel for each layer.
!                              Value must be > or = 1.
!                              UNITS:      N/A
!                              TYPE:       INTEGER
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN)
!
!       H2O_J_Index:           The "J" dimension array index position
!                              of water vapor in the input SubLevel_Absorber
!                              array - which is dimensioned Ks x J.
!                              This is necessary to properly convert the
!                              sublevel absorber amounts from ppmv to
!                              kmol.cm^-2.
!                              UNITS:      N/A
!                              TYPE:       INTEGER
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:           Character string specifying a filename
!                              in which any Messages will be logged.
!                              If not specified, or if an error occurs
!                              opening the log file, the default action
!                              is to output Messages to standard output.
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       Layer_Pressure:        Average pressure of the atmospheric layers
!                              UNITS:      hectoPascals, hPa
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Rank-1 (K-1)
!                                            K-1 == number of layers
!                              ATTRIBUTES: INTENT(IN)
!
!       Layer_Temperature:     Average temperature of the atmospheric layers
!                              UNITS:      Kelvin, K
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Same as output Layer_Pressure argument
!                              ATTRIBUTES: INTENT(IN)
!
!       Layer_Absorber:        Average absorber concentrations of the
!                              atmospheric layers
!                              UNITS:      kmol.cm^-2.
!                              TYPE:       REAL(fp)
!                              DIMENSION:  Rank-2 (K-1 x J)
!                                            K-1 == number of layers
!                                            J   == number of absorbers
!                              ATTRIBUTES: INTENT(IN)
!
! OPTIONAL OUTPUT ARGUMENTS:
!       RCS_Id:                Character string containing the Revision Control
!                              System Id field for the module.
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(OUT), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:          The return value is an integer defining the error
!                              status. The error codes are defined in the
!                              ERROR_HANDLER module.
!                              If == SUCCESS the calculation was successful
!                                 == FAILURE an unrecoverable error occurred
!                              UNITS:      N/A
!                              TYPE:       INTEGER
!                              DIMENSION:  Scalar
!
! PROCEDURE:
!       The average layer pressure is simply determined using,
!
!                                 p(k) - p(k-1)
!         Layer_Pressure(k) = --------------------
!                              LOG( p(k)/p(k-1) )
!
!       The average layer temperature is determined by summing
!       the density weighted layer temperature T.rho subLAYER
!       across the sublayers and normalising by the sum of the
!       subLAYER density,
!
!                               __ N(k)
!                              \
!                               >   Trho      [ units of kmol.cm^-2.ppmv^-1.K ]
!                              /__
!                                  1
!         Layer_Temperature = -----------
!                               __ N(k)
!                              \
!                               >   rho       [ units of kmol.cm^-2.ppmv^-1 ]
!                              /__
!                                  1
!
!
!                               __ N(k)
!                              \      1.0e-11
!                               >    --------- * dz * p
!                              /__       R
!                                  1
!                           = ---------------------------
!                               __ N(k)
!                              \       1.0e-11
!                               >     --------- * dz * p
!                              /__      R . T
!                                  1
!
!
!                               __ N(k)
!                              \
!                               >     dz . p
!                              /__
!                                  1
!                           = ----------------
!                               __ N(k)
!                              \      dz . p
!                               >    --------
!                              /__       T
!                                  1
!
!      in units of Kelvin
!
!      In the UMBC KLAYERS code, the numerator corresponds to the final
!      TSUM value (with each sublayer value corresponding to RJUNK),
!      the denominator to AJUNK, and the result to TLAY.
!
!      The average layer absorber amount is determined by simply summing
!      the sublayer absorber amount across the sublayers,
!
!                          __ N(k)
!                         \      Trho . ppmv
!        Layer_Absorber =  >    -------------
!                         /__        T
!                             1
!
!                          __ N(k)
!                         \      1.0e-11             ppmv
!                       =  >    --------- . dz . p .------
!                         /__       R                 T
!                             1
!
!       in units of kmol.cm^-2
!
!       This corresponds to ASUM (and eventually ALAY) in the 
!       UMBC KLAYERS code.
!
!       Currently, N is independent of k. That is, the same number of sublevels
!       is assumed for each layer.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 19-Jan-2001
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

  FUNCTION Integrate_Sublevels( Sublevel_Height,       &  ! Input
                                Sublevel_Pressure,     &  ! Input
                                Sublevel_Temperature,  &  ! Input
                                Sublevel_Absorber,     &  ! Input
                                n_Per_Layer,           &  ! Input
                                H2O_J_Index,           &  ! Input
                                Layer_Pressure,        &  ! Output
                                Layer_Temperature,     &  ! Output
                                Layer_Absorber,        &  ! Output
                                RCS_Id,                &  ! Revision control
                                Message_Log )          &  ! Error messaging
                              RESULT( Error_Status )
    ! Arguments
    REAL(fp),               INTENT(IN)  :: Sublevel_Height(:)
    REAL(fp),               INTENT(IN)  :: Sublevel_Pressure(:)
    REAL(fp),               INTENT(IN)  :: Sublevel_Temperature(:)
    REAL(fp),               INTENT(IN)  :: Sublevel_Absorber(:,:)
    INTEGER,                INTENT(IN)  :: n_Per_Layer
    INTEGER,                INTENT(IN)  :: H2O_J_Index
    REAL(fp),               INTENT(OUT) :: Layer_Pressure(:)
    REAL(fp),               INTENT(OUT) :: Layer_Temperature(:)
    REAL(fp),               INTENT(OUT) :: Layer_Absorber(:,:)
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: RCS_Id
    CHARACTER(*), OPTIONAL, INTENT(IN)  :: Message_Log
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Integrate_Sublevels'
    ! Local variables
    CHARACTER( 256 ) :: Message
    INTEGER :: n_Sublevels
    INTEGER :: n_Layers
    INTEGER :: n_Absorbers
    INTEGER :: j          ! Absorber index variable
    INTEGER :: k          ! Layer index variable
    INTEGER :: n, n1, n2  ! Sublevel loop/indices within *layer* k
    REAL(fp) :: Sublevel_RHOair, Sublevel_RHOair_nm1
    REAL(fp) :: Sublayer_dZ
    REAL(fp) :: Sublayer_Pressure
    REAL(fp) :: Sublayer_Temperature
    REAL(fp) :: Sublayer_T_RHOair
    REAL(fp) :: Sublayer_RHOair
    REAL(fp) :: Layer_T_RHOair_sum
    REAL(fp) :: Layer_RHOair_sum
    REAL(fp) :: Sublayer_Absorber ,Sublayer_Absorber_k
    REAL(fp) :: Layer_Absorber_sum(SIZE(Sublevel_Absorber,DIM=2))
    REAL(fp) :: Layer_dZ
    REAL(fp) :: Water_Vapor


    ! Setup
    ! -----
    Error_Status = SUCCESS
    IF ( PRESENT(RCS_Id) ) RCS_Id = MODULE_RCS_ID
    
    ! Check size of input arrays
    n_Sublevels = SIZE(Sublevel_Height)
    IF ( SIZE(Sublevel_Pressure)       < n_Sublevels .OR. &
         SIZE(Sublevel_Temperature)    < n_Sublevels .OR. &
         SIZE(Sublevel_Absorber,DIM=1) < n_Sublevels      ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Inconsistent input array sizes.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Check input array values
    IF ( ANY(Sublevel_Pressure < TOLERANCE) ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input Pressures < or = 0.0 found.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( ANY(Sublevel_Temperature < TOLERANCE) ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input Temperatures < or = 0.0 found.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    IF ( ANY(Sublevel_Absorber < ZERO) ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input absorber amounts < 0.0 found.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Check sublevel multiplier
    IF ( n_Per_Layer < 1 ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input N_PER_LAYER must be > 0.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF

    ! Check water vapour index
    IF ( H2O_J_Index < 1 .OR. &
         H2O_J_Index > SIZE(Sublevel_Absorber,DIM=2) ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Input H2O_J_Index value is invalid.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    END IF

    ! Check size of output arrays
    IF ( n_Per_Layer > 1 ) THEN
      n_Layers = n_Sublevels / n_Per_Layer
    ELSE
      n_Layers = n_Sublevels - 1
    END IF
    IF ( SIZE(Layer_Pressure)       < n_Layers .OR. &
         SIZE(Layer_Temperature)    < n_Layers .OR. &
         SIZE(Layer_Absorber,DIM=1) < n_Layers      ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Output arrays not large enough to hold result.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    n_Absorbers = SIZE(Sublevel_Absorber,DIM=2)
    IF ( SIZE(Layer_Absorber,DIM=2) < n_Absorbers ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Output Layer_Absorber array does not have '//&
                            'enough absorber dimension elements.', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF
    

    ! Calculate initial level total number density
    ! --------------------------------------------
    Sublevel_RHOair_nm1 = PP_to_ND( Sublevel_Pressure( 1 ), &
                                    Sublevel_Temperature( 1 ), &
                                    Message_Log=Message_Log )
    IF ( Sublevel_RHOair_nm1 < ZERO ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Error calculating RHOair at Sublevel 1', &
                            Error_Status, &
                            Message_Log=Message_Log )
      RETURN
    ENDIF


    ! Begin layer loop
    ! ----------------
    Layer_Loop: DO k = 1, n_Layers


      ! Initialise sum variables
      Layer_T_RHOair_sum      = ZERO
      Layer_RHOair_sum        = ZERO
      Layer_Absorber_sum( : ) = ZERO

      ! Sublevel array indices
      n1 = ( ( k-1 ) * n_Per_Layer ) + 1
      n2 = n1 + n_Per_Layer - 1


      ! Loop over sublayers
      ! -------------------
      Sublayer_Loop: DO n = n1, n2


        ! Calculate current top of sublayer total number density
        ! ------------------------------------------------------
        Sublevel_RHOair = PP_to_ND( Sublevel_Pressure( n+1 ), &
                                    Sublevel_Temperature( n+1 ), &
                                    Message_Log=Message_Log )
        IF ( Sublevel_RHOair < ZERO ) THEN
          Error_Status = FAILURE
          WRITE( Message, '( "Error calculating RHOair at Sublevel ", i4 )' ) n+1
          CALL Display_Message( ROUTINE_NAME, &
                                TRIM(Message), &
                                Error_Status, &
                                Message_Log=Message_Log )
          RETURN
        ENDIF


        ! Perform the summation for the density weighted layer temperature
        ! by summing the T.rho subLAYER product and the subLAYER density
        ! normalising factor:
        !
        !
        !                       __ N(k)
        !                      \
        !                       >   Trho      [ units of kmol.cm^-2.ppmv^-1.K ]
        !                      /__
        !                          1
        ! Layer_Temperature = -----------
        !                       __ N(k)
        !                      \
        !                       >   rho       [ units of kmol.cm^-2.ppmv^-1 ]
        !                      /__
        !                          1
        !
        !
        !                       __ N(k)
        !                      \      1.0e-11
        !                       >    --------- * dz * p
        !                      /__       R
        !                          1
        !                   = ---------------------------
        !                       __ N(k)
        !                      \       1.0e-11
        !                       >     --------- * dz * p
        !                      /__      R . T
        !                          1
        !
        !
        !                       __ N(k)
        !                      \
        !                       >     dz . p
        !                      /__
        !                          1
        !                   = ----------------
        !                       __ N(k)
        !                      \      dz . p
        !                       >    --------
        !                      /__       T
        !                          1
        !
        ! in units of Kelvin
        !
        ! In the UMBC KLAYERS code, the numerator corresponds to the final
        ! TSUM value (with each sublayer value corresponding to RJUNK),
        ! the denominator to AJUNK, and the result to TLAY.
        ! ----------------------------------------------------------------
        ! Calculate sublayer thickness, dz
        Sublayer_dZ = ABS( Sublevel_Height( n+1 ) - Sublevel_Height( n ) )

        ! Calculate sublayer pressure, p
        Sublayer_Pressure =    ( Sublevel_Pressure( n+1 ) - Sublevel_Pressure( n ) ) / &
        !                   --------------------------------------------------------
                            LOG( Sublevel_Pressure( n+1 ) / Sublevel_Pressure( n ) )

        ! Calculate sublayer temperature, T
        Sublayer_Temperature = ( Sublevel_Temperature( n+1 )*Sublevel_RHOair+Sublevel_Temperature( n )*Sublevel_RHOair_nm1 ) / &
        !                      -----------------------------------------------------------------------------------------------
                                                           ( Sublevel_RHOair+Sublevel_RHOair_nm1 )

        ! Calculate the sublayer T.rho and rho variables
        Sublayer_T_RHOair = Sublayer_dZ * Sublayer_Pressure
        Sublayer_RHOair   = Sublayer_T_RHOair / Sublayer_Temperature

        ! Sum the sublayer Trho and rho variables
        Layer_T_RHOair_sum = Layer_T_RHOair_sum + Sublayer_T_RHOair
        Layer_RHOair_sum   = Layer_RHOair_sum   + Sublayer_RHOair


        ! Perform the summation for the integrated layer absorber
        ! amount:
        !
        !                   __ N(k)
        !                  \      Trho . ppmv
        ! Layer_Absorber =  >    -------------
        !                  /__        T
        !                      1
        !
        !                   __ N(k)
        !                  \      1.0e-11             ppmv
        !                =  >    --------- . dz . p .------
        !                  /__       R                 T
        !                      1
        !
        ! in units of kmol.cm^-2
        !
        ! This corresponds to ASUM (and eventually ALAY) in the
        ! UMBC KLAYERS code.
        ! ---------------------------------------------------------
        Absorber_Sum_Loop: DO j = 1, n_Absorbers

          ! Calculate simple average sublayer absorber in ppmv
          Sublayer_Absorber = 0.5_fp * ( Sublevel_Absorber( n+1, j ) + Sublevel_Absorber( n, j ) )

          ! Convert to kmol.cm^-2
          IF ( j == H2O_J_Index ) THEN                                                  
            Sublayer_Absorber_k = PPMV_to_KMOL( Sublayer_Pressure, &
                                                Sublayer_Temperature, &
                                                Sublayer_dZ, &
                                                Sublayer_Absorber, &
                                                Message_Log=Message_Log )
            Water_Vapor = Sublayer_Absorber_k                                 
          ELSE                                                                
            Sublayer_Absorber_k = PPMV_to_KMOL( Sublayer_Pressure, &          
                                                Sublayer_Temperature, &       
                                                Sublayer_dZ, &                
                                                Sublayer_Absorber, &          
                                                Water_Vapor=Water_Vapor, &  
                                                Message_Log=Message_Log )   
          END IF                                                              

          ! Sum the column density
          Layer_Absorber_sum( j ) = Layer_Absorber_sum( j ) + Sublayer_Absorber_k

        END DO Absorber_Sum_Loop


        ! Save top boundary Density for use as bottom
        ! boundary density for next layer
        ! -------------------------------------------
        Sublevel_RHOair_nm1 = Sublevel_RHOair

      END DO Sublayer_Loop


      ! Assign the average layer values
      ! -------------------------------
      Layer_dZ = ABS( Sublevel_Height( n2+1 ) - Sublevel_Height( n1 ) )

      Layer_Pressure( k )    =    ( Sublevel_Pressure( n2+1 ) - Sublevel_Pressure( n1 ) ) / &
      !                        ----------------------------------------------------------
                               LOG( Sublevel_Pressure( n2+1 ) / Sublevel_Pressure( n1 ) )

      Layer_Temperature( k ) = Layer_T_RHOair_sum / Layer_RHOair_sum

      Layer_Absorber( k, : ) = Layer_Absorber_sum( : )

    END DO Layer_Loop

  END FUNCTION Integrate_Sublevels

END MODULE Level_Layer_Conversion
