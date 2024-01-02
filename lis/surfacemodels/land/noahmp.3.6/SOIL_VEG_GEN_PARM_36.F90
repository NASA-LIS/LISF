!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!-----------------------------------------------------------------
        SUBROUTINE SOIL_VEG_GEN_PARM_36( VEG_TBL, SOIL_TBL, GEN_TBL, MMINLU, MMINSL)
!-----------------------------------------------------------------
        USE module_sf_noahlsm_36
        use LIS_logMod, only  : LIS_logunit
        use LIS_coreMod, only : LIS_masterproc
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: VEG_TBL, SOIL_TBL, GEN_TBL, MMINLU, MMINSL
        integer :: LUMATCH, IINDEX, LC, NUM_SLOPE
        integer :: ierr
        INTEGER , PARAMETER :: OPEN_OK = 0
        character*128 :: mess , message
!-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
!             ALBBCK: SFC albedo (in percentage)
!                 Z0: Roughness length (m)
!             SHDFAC: Green vegetation fraction (in percentage)
!  Note: The ALBEDO, Z0, and SHDFAC values read from the following table
!          ALBEDO, amd Z0 are specified in LAND-USE TABLE; and SHDFAC is
!          the monthly green vegetation data
!             CMXTBL: MAX CNPY Capacity (m)
!             NROTBL: Rooting depth (layer)
!              RSMIN: Mimimum stomatal resistance (s m-1)
!              RSMAX: Max. stomatal resistance (s m-1)
!                RGL: Parameters used in radiation stress function
!                 HS: Parameter used in vapor pressure deficit functio
!               TOPT: Optimum transpiration air temperature. (K)
!             CMCMAX: Maximum canopy water capacity
!             CFACTR: Parameter used in the canopy inteception calculati
!               SNUP: Threshold snow depth (in water equivalent m) that
!                     implies 100% snow cover
!                LAI: Leaf area index (dimensionless)
!             MAXALB: Upper bound on maximum albedo over deep snow
!
!-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL
!

      OPEN(19, FILE=trim(VEG_TBL),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
      IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
        'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening VEGPARM.TBL'
        CALL wrf_error_fatal ( message )
      END IF


      LUMATCH=0

      FIND_LUTYPE : DO WHILE (LUMATCH == 0)
         READ (19,*,END=2002)
         READ (19,*,END=2002)LUTYPE
         READ (19,*)LUCATS,IINDEX

! Print messages to the lislog file instead - David Mocko
         IF(LUTYPE.EQ.MMINLU)THEN
            WRITE( mess , * ) 'LANDUSE TYPE = ' // TRIM ( LUTYPE ) // ' FOUND', LUCATS,' CATEGORIES'
!            CALL wrf_message( mess )
            if (LIS_masterproc) then
               write(LIS_logunit, *) trim(mess)
            endif
            LUMATCH=1
         ELSE
            WRITE( mess, * ) "Skipping over LUTYPE = " // TRIM ( LUTYPE )
!            call wrf_message ( "Skipping over LUTYPE = " // TRIM ( LUTYPE ) )
            if (LIS_masterproc) then
               write(LIS_logunit, *) trim(mess)
            endif
            DO LC = 1, LUCATS+12
               read(19,*)
            ENDDO
         ENDIF
      ENDDO FIND_LUTYPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
      IF ( SIZE(SHDTBL)       < LUCATS .OR. &
           SIZE(NROTBL)       < LUCATS .OR. &
           SIZE(RSTBL)        < LUCATS .OR. &
           SIZE(RGLTBL)       < LUCATS .OR. &
           SIZE(HSTBL)        < LUCATS .OR. &
           SIZE(SNUPTBL)      < LUCATS .OR. &
           SIZE(MAXALB)       < LUCATS .OR. &
           SIZE(LAIMINTBL)    < LUCATS .OR. &
           SIZE(LAIMAXTBL)    < LUCATS .OR. &
           SIZE(Z0MINTBL)     < LUCATS .OR. &
           SIZE(Z0MAXTBL)     < LUCATS .OR. &
           SIZE(ALBEDOMINTBL) < LUCATS .OR. &
           SIZE(ALBEDOMAXTBL) < LUCATS .OR. &
           SIZE(ZTOPVTBL) < LUCATS .OR. &
           SIZE(ZBOTVTBL) < LUCATS .OR. &
           SIZE(EMISSMINTBL ) < LUCATS .OR. &
           SIZE(EMISSMAXTBL ) < LUCATS ) THEN
         CALL wrf_error_fatal('Table sizes too small for value of LUCATS in module_sf_noahdrv.F')
      ENDIF

      IF(LUTYPE.EQ.MMINLU)THEN
        DO LC=1,LUCATS
            READ (19,*)IINDEX,SHDTBL(LC),                        &
                      NROTBL(LC),RSTBL(LC),RGLTBL(LC),HSTBL(LC), &
                      SNUPTBL(LC),MAXALB(LC), LAIMINTBL(LC),     &
                      LAIMAXTBL(LC),EMISSMINTBL(LC),             &
                      EMISSMAXTBL(LC), ALBEDOMINTBL(LC),         &
                      ALBEDOMAXTBL(LC), Z0MINTBL(LC), Z0MAXTBL(LC),&
    ZTOPVTBL(LC), ZBOTVTBL(LC)
        ENDDO
!
        READ (19,*)
        READ (19,*)TOPT_DATA
        READ (19,*)
        READ (19,*)CMCMAX_DATA
        READ (19,*)
        READ (19,*)CFACTR_DATA
        READ (19,*)
        READ (19,*)RSMAX_DATA
        READ (19,*)
        READ (19,*)BARE
        READ (19,*)
        READ (19,*)NATURAL
      ENDIF
!
2002   CONTINUE

      CLOSE (19)
      IF (LUMATCH == 0) then
         CALL wrf_error_fatal ("Land Use Dataset '"//MMINLU//"' not found in VEGPARM.TBL.")
      ENDIF
!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
      OPEN(19, FILE=trim(SOIL_TBL),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
      IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
        'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
        CALL wrf_error_fatal ( message )
      END IF

! Print messages to the lislog file instead - David Mocko
      WRITE(mess,*) 'INPUT SOIL TEXTURE CLASSIFICATION = ', TRIM ( MMINSL )
!      CALL wrf_message( mess )
      if (LIS_masterproc) then
         write(LIS_logunit, *) trim(mess)
      endif

      LUMATCH=0

      READ (19,*)
      READ (19,2000,END=2003)SLTYPE
2000   FORMAT (A4)
      READ (19,*)SLCATS,IINDEX
      IF(SLTYPE.EQ.MMINSL)THEN
! Print messages to the lislog file instead - David Mocko
          WRITE( mess , * ) 'SOIL TEXTURE CLASSIFICATION = ', TRIM ( SLTYPE ) , ' FOUND', &
                SLCATS,' CATEGORIES'
!          CALL wrf_message ( mess )
          if (LIS_masterproc) then
             write(LIS_logunit, *) trim(mess)
          endif
        LUMATCH=1
      ENDIF
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
      IF ( SIZE(BB    ) < SLCATS .OR. &
           SIZE(DRYSMC) < SLCATS .OR. &
           SIZE(F11   ) < SLCATS .OR. &
           SIZE(MAXSMC) < SLCATS .OR. &
           SIZE(REFSMC) < SLCATS .OR. &
           SIZE(SATPSI) < SLCATS .OR. &
           SIZE(SATDK ) < SLCATS .OR. &
           SIZE(SATDW ) < SLCATS .OR. &
           SIZE(WLTSMC) < SLCATS .OR. &
           SIZE(QTZ   ) < SLCATS  ) THEN
         CALL wrf_error_fatal('Table sizes too small for value of SLCATS in module_sf_noahdrv.F')
      ENDIF
      IF(SLTYPE.EQ.MMINSL)THEN
        DO LC=1,SLCATS
            READ (19,*) IINDEX,BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
                      REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
                      WLTSMC(LC), QTZ(LC)
        ENDDO
      ENDIF

2003   CONTINUE
      CLOSE (19)


      IF(LUMATCH.EQ.0)THEN
          CALL wrf_message( 'SOIl TEXTURE IN INPUT FILE DOES NOT ' )
          CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
          CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
      ENDIF

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
      OPEN(19, FILE=trim(GEN_TBL),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
      IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
        'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening GENPARM.TBL'
        CALL wrf_error_fatal ( message )
      END IF

      READ (19,*)
      READ (19,*)
      READ (19,*) NUM_SLOPE

        SLPCATS=NUM_SLOPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
        IF ( SIZE(slope_data) < NUM_SLOPE ) THEN
          CALL wrf_error_fatal('NUM_SLOPE too large for slope_data array in module_sf_noahdrv')
        ENDIF

        DO LC=1,SLPCATS
            READ (19,*)SLOPE_DATA(LC)
        ENDDO

        READ (19,*)
        READ (19,*)SBETA_DATA
        READ (19,*)
        READ (19,*)FXEXP_DATA
        READ (19,*)
        READ (19,*)CSOIL_DATA
        READ (19,*)
        READ (19,*)SALP_DATA
        READ (19,*)
        READ (19,*)REFDK_DATA
        READ (19,*)
        READ (19,*)REFKDT_DATA
        READ (19,*)
        READ (19,*)FRZK_DATA
        READ (19,*)
        READ (19,*)ZBOT_DATA
        READ (19,*)
        READ (19,*)CZIL_DATA
        READ (19,*)
        READ (19,*)SMLOW_DATA
        READ (19,*)
        READ (19,*)SMHIGH_DATA
        READ (19,*)
        READ (19,*)LVCOEF_DATA
      CLOSE (19)
!-----------------------------------------------------------------
      END SUBROUTINE SOIL_VEG_GEN_PARM_36 
!-----------------------------------------------------------------
