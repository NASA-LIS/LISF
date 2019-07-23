module module_sf_noahmpglacier3911

  public          :: calc_declin
  public          :: snow_init
  
  contains
    
    SUBROUTINE calc_declin( nowdate, latitude, longitude, cosz, yearlen, julian)
      use kwm_date_utilities_noahgl3911
!---------------------------------------------------------------------
      IMPLICIT NONE
!---------------------------------------------------------------------
      
      REAL, PARAMETER :: DEGRAD = 3.14159265/180.
      REAL, PARAMETER :: DPD    = 360./365.
      ! !ARGUMENTS:
      character(len=19), intent(in)  :: nowdate    ! YYYY-MM-DD_HH:mm:ss
      real,              intent(in)  :: latitude
      real,              intent(in)  :: longitude
      real,              intent(out) :: cosz
      integer,           intent(out) :: yearlen
      real,              intent(out) :: JULIAN

      REAL                           :: hrang
      real                           :: DECLIN
      real                           :: tloctim
      REAL                           :: OBECL
      REAL                           :: SINOB
      REAL                           :: SXLONG
      REAL                           :: ARG
      integer                        :: iyear
      integer                        :: iday
      integer                        :: ihour
      integer                        :: iminute
      integer                        :: isecond

      !
      ! Determine the number of days in the year
      !

      read(nowdate(1:4), '(I4)') iyear
      yearlen = 365
      if (mod(iyear,4) == 0) then
         yearlen = 366
         if (mod(iyear,100) == 0) then
            yearlen = 365
            if (mod(iyear,400) == 0) then
               yearlen = 366
               if (mod(iyear,3600) == 0) then
                  yearlen = 365
               endif
            endif
         endif
      endif

      !
      ! Determine the Julian time (floating-point day of year).
      !

      call geth_idts(nowdate(1:10), nowdate(1:4)//"-01-01", iday)
      read(nowdate(12:13), *) ihour
      read(nowdate(15:16), *) iminute
      read(nowdate(18:19), *) isecond
      julian = real(iday) + real(ihour)/24.

      !
      ! for short wave radiation

      DECLIN=0.

      !-----OBECL : OBLIQUITY = 23.5 DEGREE.

      OBECL=23.5*DEGRAD
      SINOB=SIN(OBECL)

      !-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

      IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)*DEGRAD
      IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)*DEGRAD
      ARG=SINOB*SIN(SXLONG)
      DECLIN=ASIN(ARG)

      TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + LONGITUDE/15.0 ! Local time in hours
      tloctim = AMOD(tloctim+24.0, 24.0)
      HRANG=15.*(TLOCTIM-12.)*DEGRAD
      COSZ=SIN(LATITUDE*DEGRAD)*SIN(DECLIN)+COS(LATITUDE*DEGRAD)*COS(DECLIN)*COS(HRANG)
      COSZ=MIN(COSZ,1.0);   !Added by kwH 3/1/16 to address floating point roundoff errors 
      COSZ=MAX(COSZ,-1.0);  !

      !KWM   write(wrf_err_message,10)DECDEG/DEGRAD
      !KWM10 FORMAT(1X,'*** SOLAR DECLINATION ANGLE = ',F6.2,' DEGREES.',' ***')
      !KWM   CALL wrf_debug (50, wrf_err_message)

    END SUBROUTINE calc_declin

!------------------------------------------------------------------------------------------

  SUBROUTINE SNOW_INIT ( ims , ime , jms , jme , its , itf , jts , jtf ,                  &
       &                 NSNOW , NSOIL , ZSOIL , SWE , TGXY , SNODEP ,                    &
       &                 ZSNSOXY , TSNOXY , SNICEXY ,SNLIQXY , ISNOWXY )
!------------------------------------------------------------------------------------------
!   Initialize snow arrays for Noah-MP LSM, based in input SNOWDEP, NSNOW
!   ISNOWXY is an index array, indicating the index of the top snow layer.  Valid indices
!           for snow layers range from 0 (no snow) and -1 (shallow snow) to (-NSNOW)+1 (deep snow).
!   TSNOXY holds the temperature of the snow layer.  Snow layers are initialized with 
!          temperature = ground temperature [?].  Snow-free levels in the array have value 0.0
!   SNICEXY is the frozen content of a snow layer.  Initial estimate based on SNODEP and SWE
!   SNLIQXY is the liquid content of a snow layer.  Initialized to 0.0
!   ZNSNOXY is the layer depth from the surface.  
!------------------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------------------
    INTEGER, INTENT(IN)                              :: ims, ime, jms, jme
    INTEGER, INTENT(IN)                              :: its, itf, jts, jtf
    INTEGER, INTENT(IN)                              :: NSNOW
    INTEGER, INTENT(IN)                              :: NSOIL
    REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: SWE 
    REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: SNODEP
    REAL,    INTENT(IN), DIMENSION(ims:ime, jms:jme) :: TGXY
    REAL,    INTENT(IN), DIMENSION(1:NSOIL)          :: ZSOIL

    INTEGER, INTENT(OUT), DIMENSION(ims:ime, jms:jme)                :: ISNOWXY ! Top snow layer index
    REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:NSOIL,jms:jme) :: ZSNSOXY ! Snow/soil layer depth from surface [m]
    REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:    0,jms:jme) :: TSNOXY  ! Snow layer temperature [K]
    REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:    0,jms:jme) :: SNICEXY ! Snow layer ice content [mm]
    REAL,    INTENT(OUT), DIMENSION(ims:ime, -NSNOW+1:    0,jms:jme) :: SNLIQXY ! snow layer liquid content [mm]

! Local variables:
!   DZSNO   holds the thicknesses of the various snow layers.
!   DZSNOSO holds the thicknesses of the various soil/snow layers.
    INTEGER                           :: I,J,IZ
    REAL,   DIMENSION(-NSNOW+1:    0) :: DZSNO
    REAL,   DIMENSION(-NSNOW+1:NSOIL) :: DZSNSO

!------------------------------------------------------------------------------------------

    DO J = jts , jtf
       DO I = its , itf
          IF ( SNODEP(I,J) < 0.025 ) THEN
             ISNOWXY(I,J) = 0
             DZSNO(-NSNOW+1:0) = 0.
          ELSE
             IF ( ( SNODEP(I,J) >= 0.025 ) .AND. ( SNODEP(I,J) <= 0.05 ) ) THEN
                ISNOWXY(I,J)    = -1
                DZSNO(0)  = SNODEP(I,J)
             ELSE IF ( ( SNODEP(I,J) > 0.05 ) .AND. ( SNODEP(I,J) <= 0.10 ) ) THEN
                ISNOWXY(I,J)    = -2
                DZSNO(-1) = SNODEP(I,J)/2.
                DZSNO( 0) = SNODEP(I,J)/2.
             ELSE IF ( (SNODEP(I,J) > 0.10 ) .AND. ( SNODEP(I,J) <= 0.25 ) ) THEN
                ISNOWXY(I,J)    = -2
                DZSNO(-1) = 0.05
                DZSNO( 0) = SNODEP(I,J) - DZSNO(-1)
             ELSE IF ( ( SNODEP(I,J) > 0.25 ) .AND. ( SNODEP(I,J) <= 0.45 ) ) THEN
                ISNOWXY(I,J)    = -3
                DZSNO(-2) = 0.05
                DZSNO(-1) = 0.5*(SNODEP(I,J)-DZSNO(-2))
                DZSNO( 0) = 0.5*(SNODEP(I,J)-DZSNO(-2))
             ELSE IF ( SNODEP(I,J) > 0.45 ) THEN
                ISNOWXY(I,J)     = -3
                DZSNO(-2) = 0.05
                DZSNO(-1) = 0.20
                DZSNO( 0) = SNODEP(I,J) - DZSNO(-1) - DZSNO(-2)
             ELSE
                CALL wrf_error_fatal("Problem with the logic assigning snow layers.")
             END IF
          END IF

          TSNOXY (I,-NSNOW+1:0,J) = 0.
          SNICEXY(I,-NSNOW+1:0,J) = 0.
          SNLIQXY(I,-NSNOW+1:0,J) = 0.
          DO IZ = ISNOWXY(I,J)+1 , 0
             TSNOXY(I,IZ,J)  = TGXY(I,J)  ! [k]
             SNLIQXY(I,IZ,J) = 0.00
             SNICEXY(I,IZ,J) = 1.00 * DZSNO(IZ) * (SWE(I,J)/SNODEP(I,J))  ! [kg/m3]
          END DO

          ! Assign local variable DZSNSO, the soil/snow layer thicknesses, for snow layers
          DO IZ = ISNOWXY(I,J)+1 , 0
             DZSNSO(IZ) = -DZSNO(IZ)
          END DO

          ! Assign local variable DZSNSO, the soil/snow layer thicknesses, for soil layers
          DZSNSO(1) = ZSOIL(1)
          DO IZ = 2 , NSOIL
             DZSNSO(IZ) = (ZSOIL(IZ) - ZSOIL(IZ-1))
          END DO

          ! Assign ZSNSOXY, the layer depths, for soil and snow layers
          ZSNSOXY(I,ISNOWXY(I,J)+1,J) = DZSNSO(ISNOWXY(I,J)+1)
          DO IZ = ISNOWXY(I,J)+2 , NSOIL
             ZSNSOXY(I,IZ,J) = ZSNSOXY(I,IZ-1,J) + DZSNSO(IZ)
          ENDDO

       END DO
    END DO

  END SUBROUTINE SNOW_INIT


  end module module_sf_noahmpglacier3911
