!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!-----------------------------------------------------------------
        SUBROUTINE SOIL_VEG_GEN_PARM_71( VEG_TBL, SOIL_TBL, GEN_TBL, MMINLU, MMINSL)
!-----------------------------------------------------------------
        USE module_sf_aclsm_71
        use LIS_logMod, only  : LIS_logunit
        use LIS_coreMod, only : LIS_masterproc
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: VEG_TBL, SOIL_TBL, GEN_TBL, MMINLU, MMINSL
        integer :: LUMATCH, IINDEX, LC, NUM_SLOPE
        integer :: ierr
        INTEGER , PARAMETER :: OPEN_OK = 0
        character*128 :: mess , message

      !LB later add a general parameter file and vegetation (with list of different crops in GDD for
      ! Jaemin)
!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
      OPEN(19, FILE=trim(SOIL_TBL),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
      IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
        'module_sf_aclsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
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
      IF ( SIZE(WP) < SLCATS .OR. & 
           SIZE(SAT) < SLCATS .OR. & 
           SIZE(FC) < SLCATS .OR. & 
           SIZE(INFRATE) < SLCATS .OR. & 
           SIZE(SD) < SLCATS .OR. & 
           SIZE(CL) < SLCATS .OR. & 
           SIZE(SI) < SLCATS .OR. & 
           SIZE(OC) < SLCATS ) THEN
         CALL wrf_error_fatal('Table sizes too small for value of SLCATS in AC')
      ENDIF
      IF(SLTYPE.EQ.MMINSL)THEN
        DO LC=1,SLCATS
            READ (19,*) IINDEX, wp(LC), sat(LC), fc(LC), infrate(LC), &
                        sd(LC), cl(LC), si(LC), OC(LC)
        ENDDO
      ENDIF

2003   CONTINUE
      CLOSE (19)


      IF(LUMATCH.EQ.0)THEN
          CALL wrf_message( 'SOIl TEXTURE IN INPUT FILE DOES NOT ' )
          CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
          CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
      ENDIF

!-----------------------------------------------------------------
      END SUBROUTINE SOIL_VEG_GEN_PARM_71
!-----------------------------------------------------------------
