!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! CRTM_MW_Snow_SfcOptics
!
! Module to compute the surface optical properties for SNOW surfaces at
! microwave frequencies required for determining the SNOW surface
! contribution to the radiative transfer.
!
! This module is provided to allow developers to "wrap" their existing
! codes inside the provided functions to simplify integration into
! the main CRTM_SfcOptics module.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 23-Jun-2005
!                       paul.vandelst@noaa.gov
!
!       Modified by:    Banghua Yan, 03-Oct-2007
!                       Banghua.Yan@noaa.gov
!

MODULE CRTM_MW_Snow_SfcOptics

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds,                 ONLY: fp
  USE Message_Handler,            ONLY: SUCCESS
  USE CRTM_Parameters,            ONLY: ZERO, ONE
  USE CRTM_SpcCoeff,              ONLY: SC
  USE CRTM_Surface_Define,        ONLY: CRTM_Surface_type
  USE CRTM_GeometryInfo_Define,   ONLY: CRTM_GeometryInfo_type, &
                                        CRTM_GeometryInfo_GetValue
  USE CRTM_SfcOptics_Define,      ONLY: CRTM_SfcOptics_type
  USE CRTM_SensorInfo,            ONLY: WMO_AMSUA, &
                                        WMO_AMSUB, &
                                        WMO_AMSRE, &
                                        WMO_SSMI , &
                                        WMO_MSU  , &
                                        WMO_MHS  , &
                                        WMO_SSMIS
  USE NESDIS_AMSU_SNOWEM_Module,  ONLY: NESDIS_AMSU_SNOWEM
  USE NESDIS_SSMI_SNOWEM_Module,  ONLY: NESDIS_SSMI_SnowEM
  USE NESDIS_AMSRE_SNOWEM_Module, ONLY: NESDIS_AMSRE_SNOW
  USE NESDIS_MHS_SNOWEM_Module,   ONLY: NESDIS_SNOWEM_MHS
  USE NESDIS_LandEM_Module,       ONLY: NESDIS_LandEM2
  USE NESDIS_SSMIS_SnowEM_Module, ONLY: NESDIS_SSMIS_SnowEM
  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Data types
  PUBLIC :: MWSSOVariables_type
  ! Science routines
  PUBLIC :: Compute_MW_Snow_SfcOptics
  PUBLIC :: Compute_MW_Snow_SfcOptics_TL
  PUBLIC :: Compute_MW_Snow_SfcOptics_AD


  ! -----------------
  ! Module parameters
  ! -----------------
  ! RCS Id for the module
  CHARACTER(*), PRIVATE, PARAMETER :: MODULE_RCS_ID = &
  '$Id: CRTM_MW_Snow_SfcOptics.f90 7080 2010-03-15 17:19:38Z david.groff@noaa.gov $'


  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: MWSSOVariables_type
    PRIVATE
    INTEGER :: Dummy = 0
  END TYPE MWSSOVariables_type


CONTAINS


!----------------------------------------------------------------------------------
!
! NAME:
!       Compute_MW_Snow_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at microwave
!       frequencies over a snow surface.
!
!       This function is a wrapper for third party code.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_MW_Snow_SfcOptics( Surface               , &  ! Input
!                                                 GeometryInfo          , &  ! Input
!                                                 SensorIndex           , &  ! Input
!                                                 ChannelIndex          , &  ! Output     
!                                                 SfcOptics             , &  ! Output
!                                                 MWSSOVariables        , &  ! Internal variable output
!                                                 Message_Log=Message_Log )  ! Error messaging 
!
! INPUT ARGUMENTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_Surface_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:    CRTM_GeometryInfo structure containing the 
!                        view geometry information.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_GeometryInfo_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:     Sensor index id. This is a unique index associated
!                        with a (supported) sensor used to access the
!                        shared coefficient data for a particular sensor.
!                        See the ChannelIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:    Channel index id. This is a unique index associated
!                        with a (supported) sensor channel used to access the
!                        shared coefficient data for a particular sensor's
!                        channel.
!                        See the SensorIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:     Character string specifying a filename in which any
!                        messages will be logged. If not specified, or if an
!                        error occurs opening the log file, the default action
!                        is to output messages to standard output.
!                        UNITS:      None
!                        TYPE:       CHARACTER(*)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation. On input the Angle component
!                        is assumed to contain data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_SfcOptics_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
!       MWSSOVariables:  Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_MW_Snow_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       MWSSOVariables_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output SfcOptics argument is IN OUT rather
!       than just OUT as it is assumed to contain some data upon input.
!
!
!----------------------------------------------------------------------------------

  FUNCTION Compute_MW_Snow_SfcOptics( Surface     , &  ! Input
                                      GeometryInfo, &  ! Input
                                      SensorIndex , &  ! Input
                                      ChannelIndex, &  ! Input
                                      SfcOptics   , &  ! Output
                                      MWSSOV      , &  ! Internal variable output
                                      Message_Log ) &  ! Error messaging
                                    RESULT ( Error_Status )
    ! Arguments
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics
    TYPE(MWSSOVariables_type),    INTENT(IN OUT) :: MWSSOV
    CHARACTER(*), OPTIONAL,       INTENT(IN)     :: Message_Log
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*),  PARAMETER :: ROUTINE_NAME = 'Compute_MW_Snow_SfcOptics'
    REAL(fp), PARAMETER :: MSU_SNOW_TEMPERATURE_THRESHOLD = 100.0_fp  ! K
    REAL(fp), PARAMETER :: MSU_TB_THRESHOLD               =  50.0_fp  ! K
    REAL(fp), PARAMETER :: MSU_ALPHA_C                    =   0.35_fp
    REAL(fp), PARAMETER :: MSU_EMISSIVITY_THRESHOLD       =   0.6_fp
    REAL(fp), PARAMETER :: MSU_DEFAULT_EMISSIVITY         =   0.855_fp
    REAL(fp), PARAMETER :: FREQUENCY_THRESHOLD            =  80.0_fp  ! GHz
    REAL(fp), PARAMETER :: DEFAULT_EMISSIVITY             =   0.90_fp
    REAL(fp), PARAMETER :: NOT_USED(4)                    = -99.9_fp
    INTEGER,  PARAMETER :: AMSRE_V_INDEX(6) = (/1, 3, 5, 7, 9, 11/)  ! AMSRE channels with V pol.
    INTEGER,  PARAMETER :: AMSRE_H_INDEX(6) = (/2, 4, 6, 8, 10, 12/) ! AMSRE channels with H pol.
    INTEGER,  PARAMETER :: AMSUA_INDEX(4)   = (/1, 2, 3, 15/)
    INTEGER,  PARAMETER :: SSMIS_INDEX(8)   = (/13,12,14,16,15,17,18,8/)  ! With swapped polarisations
    ! Local variables
    INTEGER :: i
    REAL(fp) :: Sensor_Zenith_Angle
    REAL(fp) :: Alpha


    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS
    CALL CRTM_GeometryInfo_GetValue( GeometryInfo, Sensor_Zenith_Angle = Sensor_Zenith_Angle )


    ! --------------------------------
    ! Compute the surface emissivities
    ! --------------------------------
    Sensor_Type: SELECT CASE( Surface%SensorData%WMO_Sensor_ID )

      ! AMSU-A emissivity model
      CASE( WMO_AMSUA )                                                                                 
        DO i = 1, SfcOptics%n_Angles                                                                    
          CALL NESDIS_AMSU_SNOWEM( Sensor_Zenith_Angle,                     &  ! Input, Degree           
                                   SfcOptics%Angle(i),                      &  ! Input, Degree           
                                   SC(SensorIndex)%Frequency(ChannelIndex), &  ! Input, GHz                  
                                   Surface%Snow_Depth,                      &  ! Input, mm               
                                   Surface%Snow_Temperature,                &  ! Input, K                
                                   Surface%SensorData%Tb(AMSUA_INDEX),      &  ! Input, AMSUA            
                                   NOT_USED(1:2),                           &  ! Input, AMSUB  *** NO AMSU-B DATA ***            
                                   SfcOptics%Emissivity(i,2),               &  ! Output, H component      
                                   SfcOptics%Emissivity(i,1)                )  ! Output, V component     
        END DO                                                                                           

      ! AMSU-B emissivity model
      CASE( WMO_AMSUB)
        DO i = 1, SfcOptics%n_Angles                                                                    
          CALL NESDIS_AMSU_SNOWEM( Sensor_Zenith_Angle,                     &  ! Input, Degree           
                                   SfcOptics%Angle(i),                      &  ! Input, Degree           
                                   SC(SensorIndex)%Frequency(ChannelIndex), &  ! Input, GHz                  
                                   Surface%Snow_Depth,                      &  ! Input, mm               
                                   Surface%Snow_Temperature,                &  ! Input, K                
                                   NOT_USED,                                &  ! Input  AMSUA  *** NO AMSU-A DATA ***            
                                   Surface%SensorData%Tb(1:2),              &  ! Input, AMSUB            
                                   SfcOptics%Emissivity(i,2),               &  ! Output, H component      
                                   SfcOptics%Emissivity(i,1)                )  ! Output, V component     
        END DO                                                                                           

      ! MHS emissivity model
      CASE (WMO_MHS)
        DO i = 1, SfcOptics%n_Angles
          CALL NESDIS_SNOWEM_MHS( Sensor_Zenith_Angle,                     &  ! Input, Degree
                                  SfcOptics%Angle(i),                      &  ! Input, Degree
                                  SC(SensorIndex)%Frequency(ChannelIndex), &  ! Input, GHz
                                  Surface%Snow_Temperature,                &  ! Input, K
                                  Surface%SensorData%Tb(1:2),              &  ! Input, AMSUB
                                  SfcOptics%Emissivity(i,2),               &  ! Output, H component
                                  SfcOptics%Emissivity(i,1)                )  ! Output, V component
        END DO

      ! AMSR-E emissivity model
      CASE( WMO_AMSRE )                                                                                 
        DO i = 1, SfcOptics%n_Angles                                                                    
          CALL NESDIS_AMSRE_SNOW(SC(SensorIndex)%Frequency(ChannelIndex), &  ! Input, GHz                  
                                 SfcOptics%Angle(i),                      &  ! Input, Degree            
                                 Surface%SensorData%Tb(AMSRE_V_INDEX),    &  ! Input, Tb_V, K           
                                 Surface%SensorData%Tb(AMSRE_H_INDEX),    &  ! Input, Tb_H, K           
                                 Surface%Snow_Temperature,                &  ! Input, Ts, K             
                                 Surface%Snow_Temperature,                &  ! Input, Tsnow, K          
                                 SfcOptics%Emissivity(i,2),               &  ! Output, H component       
                                 SfcOptics%Emissivity(i,1)                )  ! Output, V component      
        END DO                                                                                           

      ! SSM/I emissivity model
      CASE( WMO_SSMI )                                                                                 
        DO i = 1, SfcOptics%n_Angles                                                                    
          CALL NESDIS_SSMI_SnowEM(SC(SensorIndex)%Frequency(ChannelIndex), &  ! Input, GHz                  
                                  SfcOptics%Angle(i),                      &  ! Input, Degree           
                                  Surface%Snow_Temperature,                &  ! Input, K                
                                  Surface%SensorData%Tb,                   &  ! Input, K                
                                  Surface%Snow_Depth,                      &  ! Input, mm               
                                  SfcOptics%Emissivity(i,2),               &  ! Output, H component      
                                  SfcOptics%Emissivity(i,1)                )  ! Output, V component     
        END DO                                                                                           

      ! SSMIS emissivity model
      CASE( WMO_SSMIS )                                                                                 
        DO i = 1, SfcOptics%n_Angles                                                                    
          CALL NESDIS_SSMIS_SnowEM(SC(SensorIndex)%Frequency(ChannelIndex), &  ! Input, GHz                  
                                   SfcOptics%Angle(i),                      &  ! Input, Degree           
                                   Surface%Snow_Temperature,                &  ! Input, K                
                                   Surface%SensorData%Tb(SSMIS_INDEX),      &  ! Input, K                
                                   Surface%Snow_Depth,                      &  ! Input, mm               
                                   SfcOptics%Emissivity(i,2),               &  ! Output, H component      
                                   SfcOptics%Emissivity(i,1)                )  ! Output, V component     
        END DO                                                                                           

      ! MSU emissivity model
      CASE( WMO_MSU )                                                                                 
        DO i = 1, SfcOptics%n_Angles  
          IF( Surface%Snow_Temperature > MSU_SNOW_TEMPERATURE_THRESHOLD .AND. &
              Surface%SensorData%Tb(1) > MSU_TB_THRESHOLD                     ) THEN
            Alpha = MSU_ALPHA_C * Surface%Snow_Temperature
            SfcOptics%Emissivity(i,1) = (Surface%SensorData%Tb(1)-Alpha)/&
                                        (Surface%Snow_Temperature-Alpha)
            IF( SfcOptics%Emissivity(i,1) > ONE ) &
              SfcOptics%Emissivity(i,1) = ONE
            IF( SfcOptics%Emissivity(i,1) < MSU_EMISSIVITY_THRESHOLD ) &
              SfcOptics%Emissivity(i,1) = MSU_EMISSIVITY_THRESHOLD
          ELSE
            SfcOptics%Emissivity(i,1) = MSU_DEFAULT_EMISSIVITY
          END IF
          SfcOptics%Emissivity(i,2) = SfcOptics%Emissivity(i,1)
        END DO

      ! Default physical model
      CASE DEFAULT
        IF ( SC(SensorIndex)%Frequency(ChannelIndex) < FREQUENCY_THRESHOLD ) THEN
          DO i = 1, SfcOptics%n_Angles
!YDT            CALL NESDIS_LandEM( SfcOptics%Angle(i),                      & ! Input, Degree
!                       SC(SensorIndex)%Frequency(ChannelIndex), & ! Input, GHz
!                      NOT_USED(1),                             & ! Input, Soil_Moisture_Content, g.cm^-3
!                                NOT_USED(1),                             & ! Input, Vegetation_Fraction
!                                Surface%Snow_Temperature,                & ! Input, K
!                                Surface%Snow_Temperature,                & ! Input, K
!                                Surface%Snow_Depth,                      & ! Input, mm
!                                SfcOptics%Emissivity(i,2),               & ! Output, H component
!                                SfcOptics%Emissivity(i,1)                ) ! Output, V component

            CALL NESDIS_LandEM2( SfcOptics%Angle(i),                      & ! Input, Degree
                                SC(SensorIndex)%Frequency(ChannelIndex), & ! Input, GHz
                                Surface,                                 & ! Input
                                SfcOptics%Emissivity(i,2),               & ! Output, H component
                                SfcOptics%Emissivity(i,1)                ) ! Output, V component

          END DO
        ELSE
          SfcOptics%Emissivity(1:SfcOptics%n_Angles,1:2) = DEFAULT_EMISSIVITY
        END IF

    END SELECT Sensor_Type

    ! -----------------------------------
    ! Compute the surface reflectivities, 
    ! assuming a specular surface
    ! -----------------------------------
    SfcOptics%Reflectivity = ZERO
    DO i = 1, SfcOptics%n_Angles
      SfcOptics%Reflectivity(i,1,i,1) = ONE-SfcOptics%Emissivity(i,1)
      SfcOptics%Reflectivity(i,2,i,2) = ONE-SfcOptics%Emissivity(i,2)
    END DO

  END FUNCTION Compute_MW_Snow_SfcOptics


!----------------------------------------------------------------------------------
!
! NAME:
!       Compute_MW_Snow_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at microwave frequencies over a snow surface.
!
!       This function is a wrapper for third party code.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_MW_Snow_SfcOptics_TL( Surface               , &  ! Input
!                                                    SfcOptics             , &  ! Input     
!                                                    Surface_TL            , &  ! Input
!                                                    GeometryInfo          , &  ! Input
!                                                    SensorIndex           , &  ! Input
!                                                    ChannelIndex          , &  ! Output     
!                                                    SfcOptics_TL          , &  ! Output
!                                                    MWSSOVariables        , &  ! Internal variable input
!                                                    Message_Log=Message_Log )  ! Error messaging 
!
! INPUT ARGUMENTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_Surface_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Surface_TL:      CRTM_Surface structure containing the tangent-linear 
!                        surface state data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_Surface_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_SfcOptics_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:    CRTM_GeometryInfo structure containing the 
!                        view geometry information.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_GeometryInfo_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:     Sensor index id. This is a unique index associated
!                        with a (supported) sensor used to access the
!                        shared coefficient data for a particular sensor.
!                        See the ChannelIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:    Channel index id. This is a unique index associated
!                        with a (supported) sensor channel used to access the
!                        shared coefficient data for a particular sensor's
!                        channel.
!                        See the SensorIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       MWSSOVariables:  Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_MW_Snow_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       MWSSOVariables_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:     Character string specifying a filename in which any
!                        messages will be logged. If not specified, or if an
!                        error occurs opening the log file, the default action
!                        is to output messages to standard output.
!                        UNITS:      None
!                        TYPE:       CHARACTER(*)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       SfcOptics_TL:    CRTM_SfcOptics structure containing the tangent-linear
!                        surface optical properties required for the tangent-
!                        linear radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_SfcOptics_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output SfcOptics_TL argument is IN OUT rather
!       than just OUT. This is necessary because the argument may be defined
!       upon input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!
!----------------------------------------------------------------------------------


  FUNCTION Compute_MW_Snow_SfcOptics_TL( Surface     , &  ! Input
                                         SfcOptics   , &  ! Input     
                                         Surface_TL  , &  ! Input
                                         GeometryInfo, &  ! Input
                                         SensorIndex , &  ! Input
                                         ChannelIndex, &  ! Input
                                         SfcOptics_TL, &  ! Output     
                                         MWSSOV      , &  ! Internal variable input
                                         Message_Log ) &  ! Error messaging 
                                       RESULT ( Error_Status )
    ! Arguments
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface_TL
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_TL
    TYPE(MWSSOVariables_type),    INTENT(IN)     :: MWSSOV
    CHARACTER(*), OPTIONAL,       INTENT(IN)     :: Message_Log
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Snow_SfcOptics_TL'
    ! Local variables


    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS


    ! -----------------------------------------------------
    ! Compute the tangent-linear surface optical parameters
    !
    ! ***No TL models yet, so default TL output is zero***
    ! -----------------------------------------------------
    SfcOptics_TL%Reflectivity = ZERO
    SfcOptics_TL%Emissivity   = ZERO

  END FUNCTION Compute_MW_Snow_SfcOptics_TL



!----------------------------------------------------------------------------------
!
! NAME:
!       Compute_MW_Snow_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at microwave frequencies over a snow surface.
!
!       This function is a wrapper for third party code.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_MW_Snow_SfcOptics_AD( Surface               , &  ! Input
!                                                    SfcOptics             , &  ! Input     
!                                                    SfcOptics_AD          , &  ! Input     
!                                                    GeometryInfo          , &  ! Input
!                                                    SensorIndex           , &  ! Input
!                                                    ChannelIndex          , &  ! Output     
!                                                    Surface_AD            , &  ! Output
!                                                    MWSSOVariables        , &  ! Internal variable input
!                                                    Message_Log=Message_Log )  ! Error messaging 
!
! INPUT ARGUMENTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_Surface_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_SfcOptics_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SfcOptics_AD:    CRTM_SfcOptics structure containing the adjoint
!                        surface optical properties required for the adjoint
!                        radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_SfcOptics_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
!       GeometryInfo:    CRTM_GeometryInfo structure containing the 
!                        view geometry information.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_GeometryInfo_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:     Sensor index id. This is a unique index associated
!                        with a (supported) sensor used to access the
!                        shared coefficient data for a particular sensor.
!                        See the ChannelIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:    Channel index id. This is a unique index associated
!                        with a (supported) sensor channel used to access the
!                        shared coefficient data for a particular sensor's
!                        channel.
!                        See the SensorIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       MWSSOVariables:  Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_MW_Snow_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       MWSSOVariables_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Message_Log:     Character string specifying a filename in which any
!                        messages will be logged. If not specified, or if an
!                        error occurs opening the log file, the default action
!                        is to output messages to standard output.
!                        UNITS:      None
!                        TYPE:       CHARACTER(*)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       Surface_AD:      CRTM_Surface structure containing the adjoint
!                        surface state data.
!                        UNITS:      N/A
!                        TYPE:       TYPE(CRTM_Surface_type)
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the input SfcOptics_AD argument is IN OUT rather
!       than just OUT. This is necessary because components of this argument
!       may need to be zeroed out upon output.
!
!       Note the INTENT on the output Surface_AD argument is IN OUT rather
!       than just OUT. This is necessary because the argument may be defined
!       upon input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!
!----------------------------------------------------------------------------------

  FUNCTION Compute_MW_Snow_SfcOptics_AD( Surface     , &  ! Input
                                         SfcOptics   , &  ! Input     
                                         SfcOptics_AD, &  ! Input
                                         GeometryInfo, &  ! Input
                                         SensorIndex , &  ! Input
                                         ChannelIndex, &  ! Input
                                         Surface_AD  , &  ! Output     
                                         MWSSOV      , &  ! Internal variable input
                                         Message_Log ) &  ! Error messaging 
                                       RESULT ( Error_Status )
    ! Arguments
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_AD
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_Surface_type),      INTENT(IN OUT) :: Surface_AD
    TYPE(MWSSOVariables_type),    INTENT(IN)     :: MWSSOV
    CHARACTER(*), OPTIONAL,       INTENT(IN)     :: Message_Log
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Snow_SfcOptics_AD'
    ! Local variables


    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS


    ! ----------------------------------------------
    ! Compute the adjoint surface optical parameters
    !
    ! ***No AD models yet, so there is no impact on AD result***
    ! ----------------------------------------------
    SfcOptics_AD%Reflectivity = ZERO
    SfcOptics_AD%Emissivity   = ZERO

  END FUNCTION Compute_MW_Snow_SfcOptics_AD

END MODULE CRTM_MW_Snow_SfcOptics
