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
SUBROUTINE SOIL_PARM_AC72(SOIL_TBL)
  !-----------------------------------------------------------------
  use LIS_logMod, only  : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_endrun
  use LIS_coreMod, only : LIS_masterproc
  USE module_sf_aclsm_72
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: SOIL_TBL
  integer :: IINDEX, LC
  integer :: ierr
  INTEGER , PARAMETER :: OPEN_OK = 0
  integer :: ftn
  character*128 :: mess , message

  !-----SPECIFY SYSTEM PARAMETERS:
  !
  !-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
  ftn = LIS_getNextUnitNumber()
  OPEN(ftn, FILE=trim(SOIL_TBL),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'SOIL_PARM_AC72.F90: soil_veg_gen_parm: failure opening SOILPARM.TBL'
     write(LIS_logunit,*) trim(message)
     call LIS_endrun
  END IF

  READ (ftn,*)
  READ (ftn,*)SLCATS,IINDEX
  WRITE( mess , * ) 'AC72: ', SLCATS,' CATEGORIES'
  if (LIS_masterproc) then
     write(LIS_logunit, *) trim(mess)
  endif
  ! prevent possible array overwrite
  IF ( SIZE(WP) < SLCATS .OR. & 
       SIZE(SAT) < SLCATS .OR. & 
       SIZE(FC) < SLCATS .OR. & 
       SIZE(INFRATE) < SLCATS .OR. & 
       SIZE(SD) < SLCATS .OR. & 
       SIZE(CL) < SLCATS .OR. & 
       SIZE(SI) < SLCATS .OR. & 
       SIZE(OC) < SLCATS ) THEN
     write(LIS_logunit,*) 'Table sizes too small for value of SLCATS in AC72'
     call LIS_endrun
  ENDIF
  DO LC=1,SLCATS
     READ (ftn,*) IINDEX, wp(LC), sat(LC), fc(LC), infrate(LC), &
          sd(LC), cl(LC), si(LC), OC(LC)
  ENDDO
  call LIS_releaseUnitNumber(ftn)

  !-----------------------------------------------------------------
END SUBROUTINE SOIL_PARM_AC72
!-----------------------------------------------------------------
