!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!-----------------------------------------------------------------
        SUBROUTINE SOIL_PARM_72(SOIL_TBL)
!-----------------------------------------------------------------
        USE module_sf_aclsm_72
        use LIS_logMod, only  : LIS_logunit
        use LIS_coreMod, only : LIS_masterproc
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: SOIL_TBL
        integer :: IINDEX, LC
        integer :: ierr
        INTEGER , PARAMETER :: OPEN_OK = 0
        character*128 :: mess , message

!-----SPECIFY SYSTEM PARAMETERS:
!LB: other parameters could be added  in the future
!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
      OPEN(19, FILE=trim(SOIL_TBL),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
      IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
        'module_sf_aclsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
        CALL wrf_error_fatal ( message )
      END IF

      READ (19,*)
      READ (19,*)SLCATS,IINDEX
      WRITE( mess , * ) 'AC72: ', SLCATS,' CATEGORIES'
      if (LIS_masterproc) then
            write(LIS_logunit, *) trim(mess)
      endif
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
      IF ( SIZE(WP) < SLCATS .OR. & 
           SIZE(SAT) < SLCATS .OR. & 
           SIZE(FC) < SLCATS .OR. & 
           SIZE(INFRATE) < SLCATS .OR. & 
           SIZE(SD) < SLCATS .OR. & 
           SIZE(CL) < SLCATS .OR. & 
           SIZE(SI) < SLCATS .OR. & 
           SIZE(OC) < SLCATS ) THEN
         CALL wrf_error_fatal('Table sizes too small for value of SLCATS in AC72')
      ENDIF
      DO LC=1,SLCATS
      READ (19,*) IINDEX, wp(LC), sat(LC), fc(LC), infrate(LC), &
                  sd(LC), cl(LC), si(LC), OC(LC)
      ENDDO
      CLOSE (19)

!-----------------------------------------------------------------
      END SUBROUTINE SOIL_PARM_72
!-----------------------------------------------------------------
