!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT)
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! SUBROUTINE: ARFSSMRETRIEVAL
!
! REVISION HISTORY:
!  22 Feb 2022: P.W.LIU; Initial implemetation
!  22 Feb 2022: Yonghwan Kwon; modified for LDT
!
! DESCRIPTION: RETRIEVE SMAP SM FOR ARFS
! INPUT : SMAP - L1B Brightness Temperature               
! OUTPUT: SMAPTB_ARFSGRIDE_ddmmyyy.dat
! NOTES : Inverse Distance Squared with 0.4 deg serching window
!-------------------------------------------------------------------------

 subroutine ARFSSMRETRIEVAL(SMAPFILE,TS_bfresample_01,TS_bfresample_02,&
                            ARFS_SNOW,DOY,UTChr)

    !USE HDF5
    USE VARIABLES
    USE DATADOMAIN
    USE FUNCTIONS
    USE netcdf
    USE invdist_temp2smap
    USE varsio_m
    USE algo_vpol_m
    use LDT_logMod         
    USE LDT_smap_e_oplMod  
 
    IMPLICIT NONE
! !ARGUMENTS:
    CHARACTER (len=100)          :: SMAPFILE                             
    REAL*4, DIMENSION(2560,1920) :: TS_bfresample_01, TS_bfresample_02   
    REAL*4, DIMENSION(2560,1920) :: ARFS_SNOW, UTChr                     
    INTEGER                      :: DOY                                  
!EOP 
    INTEGER :: i, j, nrow, mcol          
    CHARACTER (len=100) :: fname_TAU    
    CHARACTER (len=5) :: DOY_chr
    REAL*4 :: C, K, sm_retrieval, tau_return
    REAL*4, DIMENSION(2560,1920) :: ARFS_TB
    REAL*4, DIMENSION(2560,1920) :: ARFS_TAU, ARFS_CLAY, ARFS_BD, ARFS_OMEGA, ARFS_H
    INTEGER*1, DIMENSION(2560,1920) :: ARFS_LC, ARFS_SM_FLAG
    INTEGER*1 :: retrieval_flag
    REAL*4, DIMENSION(2560,1920) :: ARFS_TS_01, ARFS_TS_02, ARFS_SM
    REAL*8 ,DIMENSION(:), ALLOCATABLE :: ARFS_FINE_LAT, ARFS_FINE_LON
    REAL*8 ,DIMENSION(:), ALLOCATABLE :: ARFS_LAT, ARFS_LON
    REAL :: T1, T2
    
    INTEGER*4 :: ios, NX, NY
    INTEGER*4 :: ncid, nid, tsoil01id

    character (len=100) :: retrieval_fname
    integer             :: L1B_dir_len,L1B_fname_len
    real                :: utc_check

    nrow=2560
    mcol=1920

    ! RESAMPLE TEFF TO 33 KM on ARFS GRID
    write (LDT_logunit,*) '[INFO] Resampling effective soil temperature'
    CALL ARFS_GEO
    ALLOCATE(ARFS_LAT(arfs_nrow_lat),ARFS_LON(arfs_mcol_lon))
    ARFS_LAT = LAT(arfs_geo_lat_lo,arfs_geo_lat_up,-arfs_lat_space)
    ARFS_LON = LON(arfs_geo_lon_lf,arfs_geo_lon_rt,arfs_lon_space)
    CALL ARFS_3KM_GEO
    ALLOCATE(ARFS_FINE_LAT(arfs_nrow_lat_3km),ARFS_FINE_LON(arfs_mcol_lon_3km))
    ARFS_FINE_LAT = LAT(arfs_geo_lat_lo,arfs_geo_lat_up,-arfs_lat_3km_space)
    ARFS_FINE_LON = LON(arfs_geo_lon_lf,arfs_geo_lon_rt,arfs_lon_3km_space)
    CALL RESAMPLETEMP(TS_bfresample_01,ARFS_LAT,ARFS_LON,ARFS_FINE_LAT,ARFS_FINE_LON,ARFS_TS_01)
    CALL RESAMPLETEMP(TS_bfresample_02,ARFS_LAT,ARFS_LON,ARFS_FINE_LAT,ARFS_FINE_LON,ARFS_TS_02)
    ! IF EVENTUALLY THE RESAMPLING DOES NOT CHANGE TEFF MUCH WE COULD SIMPLELY USE ARFS_TS=TS_bfresample
    ! UP TO HERE TAKES 38 SECS
    write (LDT_logunit,*) '[INFO] Finished resampling effective soil temperature'

    ! get RESAMPLED TB
    ARFS_TB = SMAPeOPL%ARFS_TBV_COR

    ! LOAD TAU ------------------------------------------------------------
    write(DOY_chr,"(I0.3)") DOY
    fname_TAU = trim(SMAPeOPL%TAUdir)//"/tau_cmg_arfs_"//trim(DOY_chr)//".dat"
    write (LDT_logunit,*) '[INFO] Reading TAU from ', trim(fname_TAU)
    OPEN(UNIT=1,FILE=fname_TAU,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*nrow*mcol,STATUS='OLD',convert='little_endian')
    READ(1, rec=1) ARFS_TAU
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading TAU'

    ! LOAD OMEGA
    write (LDT_logunit,*) '[INFO] Reading OMEGA from ', trim(SMAPeOPL%OMEGAfile)
    OPEN(UNIT=1,FILE=SMAPeOPL%OMEGAfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*nrow*mcol,STATUS='OLD', convert='little_endian')
    READ(1, rec=1) ARFS_OMEGA
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading OMEGA'

    ! LOAD SOIL
    write (LDT_logunit,*) '[INFO] Reading soil bulk density from ', trim(SMAPeOPL%BDfile)
    OPEN(UNIT=1,FILE=SMAPeOPL%BDfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*nrow*mcol,STATUS='OLD', convert='little_endian') !Bulk Density
    READ(1, rec=1) ARFS_BD
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading soil bulk density'

    write (LDT_logunit,*) '[INFO] Reading soil clay fraction from ', trim(SMAPeOPL%CLAYfile)
    OPEN(UNIT=1,FILE=SMAPeOPL%CLAYfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*nrow*mcol,STATUS='OLD',convert='little_endian') !Clay Fraction
    READ(1, rec=1) ARFS_CLAY
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading soil clay fraction'

    ! LOAD ROUGHNESS
    write (LDT_logunit,*) '[INFO] Reading roughness from ', trim(SMAPeOPL%Hfile)
    OPEN(UNIT=1,FILE=SMAPeOPL%Hfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*nrow*mcol,STATUS='OLD',convert='little_endian') !roughness
    READ(1, rec=1) ARFS_H
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading roughness'

    ! LOAD LANDCOVER
    write (LDT_logunit,*) '[INFO] Reading landcover from ', trim(SMAPeOPL%LCfile)
    OPEN(UNIT=1,FILE=SMAPeOPL%LCfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1*nrow*mcol,STATUS='OLD',convert='little_endian')
    READ(1, rec=1) ARFS_LC
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading landcover'

    !generate soil moisture retrievals
    write (LDT_logunit,*) '[INFO] Generating soil moisture retrievals' 
    ARFS_SM=-9999
    ARFS_SM_FLAG=-1
    DO i=1,nrow !ROW LON
       DO j=1,mcol !COL LAT
          tbv = ARFS_TB(i,j)

          if(UTChr(i,j).ge.0) then
             utc_check = UTChr(i,j) - floor(UTChr(i,j))

             if(utc_check.le.0.5) then
                Ts  = ARFS_TS_01(i,j)
             else
                Ts  = ARFS_TS_02(i,j)
             endif
          endif

          IF (tbv.GT.0.0.AND.Ts.GT.0.AND.ARFS_SNOW(i,j).LE.SMAPeOPL%SD_thold.AND.ARFS_BD(i,j).NE.-9999.AND.ARFS_LC(i,j).NE.0.AND.&
            UTChr(i,j).GE.0) THEN
             bulkdensity = ARFS_BD(i,j)
             clay = ARFS_CLAY(i,j)
             tau = ARFS_TAU(i,j)
             omega = ARFS_OMEGA(i,j)
             h = ARFS_H(i,j)
             topigbptype = ARFS_LC(i,j)

             CALL algo_vpol(real(i),real(j),sm_retrieval, tau_return, retrieval_flag)
             ARFS_SM(i,j)=sm_retrieval
             ARFS_SM_FLAG(i,j)=retrieval_flag
          ELSE
             !PRINT*,i, j, "NO RETRIEVAL"
          END IF
       END DO !jj=1,mcol !COL LAT
    END DO !ii=1,nrow !ROW LON

    !write soil moisture retrieval outputs
    L1B_dir_len = len_trim(SMAPeOPL%L1Bdir)
    L1B_fname_len = len_trim(SMAPFILE)

    if(SMAPeOPL%L1Btype.eq.1) then  !NRT
       retrieval_fname = trim(SMAPeOPL%SMoutdir)//"/"//"ARFS_SM_V_"//&
                         trim(SMAPFILE(L1B_dir_len+18:L1B_fname_len-3))//".dat"
    elseif(SMAPeOPL%L1Btype.eq.2) then  !Historical
       retrieval_fname = trim(SMAPeOPL%SMoutdir)//"/"//"ARFS_SM_V_"//&
                         trim(SMAPFILE(L1B_dir_len+14:L1B_fname_len-3))//".dat"
    endif

    write (LDT_logunit,*) '[INFO] Writing soil moisture retrieval file ', trim(retrieval_fname)
    OPEN(UNIT=151, FILE=retrieval_fname,FORM='UNFORMATTED',ACCESS='DIRECT', RECL=nrow*mcol*4)
    WRITE(UNIT=151, REC = 1) ARFS_SM
    CLOSE(151)
    write (LDT_logunit,*) '[INFO] Successfully wrote soil moisture retrieval file ', trim(retrieval_fname)
    write (LDT_logunit,*) '[INFO] Finished generating soil moisture retrievals'

 end subroutine ARFSSMRETRIEVAL
