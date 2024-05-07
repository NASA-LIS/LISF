!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! SUBROUTINE: ARFSSMRETRIEVAL
!
! REVISION HISTORY:
!  22 Feb 2022: P.W.LIU; Initial implemetation
!  22 Feb 2022: Yonghwan Kwon; modified for LDT
!  10 Feb 2023: Eric Kemp, modified to output retrievals in netCDF.
!  21 Feb 2023: Eric Kemp, added third LIS time level.
!
! DESCRIPTION: RETRIEVE SMAP SM FOR ARFS
! INPUT : SMAP - L1B Brightness Temperature
! OUTPUT: SMAPTB_ARFSGRIDE_ddmmyyy.dat
! NOTES : Inverse Distance Squared with 0.4 deg serching window
!-------------------------------------------------------------------------

subroutine ARFSSMRETRIEVAL(SMAPFILE, &
     TS_bfresample_01, TS_bfresample_02, TS_bfresample_03, &
     ARFS_SNOW, DOY, UTChr, firsttime, secondtime, thirdtime)

   !USE HDF5
    use esmf
    USE VARIABLES
    USE DATADOMAIN
    USE FUNCTIONS
    USE netcdf
    USE invdist_temp2smap
    USE varsio_m
    USE algo_vpol_m
    use LDT_ARFSSM_netcdfMod, only: LDT_ARFSSM_write_netcdf
    use LDT_logMod, only: LDT_logunit
    USE LDT_smap_e_oplMod  

    IMPLICIT NONE
! !ARGUMENTS:
    CHARACTER (len=100)          :: SMAPFILE                             
    REAL*4, DIMENSION(2560,1920), intent(in) :: TS_bfresample_01, &
         TS_bfresample_02, TS_bfresample_03
    REAL*4, DIMENSION(2560,1920) :: ARFS_SNOW, UTChr                     
    INTEGER                      :: DOY                                  
    type(ESMF_Time), intent(in) :: firsttime
    type(ESMF_Time), intent(in) :: secondtime
    type(ESMF_Time), intent(in) :: thirdtime
!EOP 
    INTEGER :: i, j, nrow, mcol          
    CHARACTER (len=100) :: fname_TAU    
    CHARACTER (len=5) :: DOY_chr
    REAL*4 :: C, K, sm_retrieval, tau_return
    REAL*4, DIMENSION(2560,1920) :: ARFS_TB
    REAL*4, DIMENSION(2560,1920) :: ARFS_TAU, ARFS_CLAY, ARFS_BD, ARFS_OMEGA, ARFS_H
    INTEGER*1, DIMENSION(2560,1920) :: ARFS_LC, ARFS_SM_FLAG
    INTEGER*1 :: retrieval_flag
    REAL*4, DIMENSION(2560,1920) :: ARFS_TS_01, ARFS_TS_02, ARFS_TS_03, ARFS_SM
    REAL*8 ,DIMENSION(:), ALLOCATABLE :: ARFS_FINE_LAT, ARFS_FINE_LON
    REAL*8 ,DIMENSION(:), ALLOCATABLE :: ARFS_LAT, ARFS_LON
    REAL :: T1, T2
    
    INTEGER*4 :: ios, NX, NY
    INTEGER*4 :: ncid, nid, tsoil01id

    character (len=100) :: retrieval_fname
    integer             :: L1B_dir_len,L1B_fname_len
    real                :: utc_check

    ! EMK
    character(8) :: yyyymmdd
    character(6) :: hhmmss
    real :: deltasec, wgt
    integer :: firstUTCyr, firstUTCmo, firstUTCdy, firstUTChr
    integer :: secondUTCyr, secondUTCmo, secondUTCdy, secondUTChr
    integer :: thirdUTCyr, thirdUTCmo, thirdUTCdy, thirdUTChr

    real :: TS_A, TS_B

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
    CALL RESAMPLETEMP(TS_bfresample_03,ARFS_LAT,ARFS_LON,ARFS_FINE_LAT,ARFS_FINE_LON,ARFS_TS_03)
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

    call ESMF_TimeGet(firsttime, yy=firstUTCyr, mm=firstUTCmo, dd=firstUTCdy, &
         h=firstUTChr)
    call ESMF_TimeGet(secondtime, yy=secondUTCyr, mm=secondUTCmo, &
         dd=secondUTCdy, &
         h=secondUTChr)
    call ESMF_TimeGet(thirdtime, yy=thirdUTCyr, mm=thirdUTCmo, dd=thirdUTCdy, &
         h=thirdUTChr)

    DO j=1,mcol !COL LAT
       DO i=1,nrow !ROW LON

          tbv = ARFS_TB(i,j)

          if (UTChr(i,j) < 0) cycle

          if (UTChr(i,j) == firstUTChr) then
             TS_A = ARFS_TS_01(i,j)
             TS_B = ARFS_TS_02(i,j)
             wgt = 1
          else if (UTChr(i,j) > firstUTChr .and. &
               UTChr(i,j) < secondUTChr) then
             TS_A = ARFS_TS_01(i,j)
             TS_B = ARFS_TS_02(i,j)
             deltasec = ( UTChr(i,j) - firstUTChr ) * 3600
             wgt = (10800. - deltasec) / 10800.
          else if (UTChr(i,j) > firstUTChr .and. &
               firstUTChr == 21 .and. secondUTChr == 0) then
             TS_A = ARFS_TS_01(i,j)
             TS_B = ARFS_TS_02(i,j)
             deltasec = ( UTChr(i,j) - firstUTChr ) * 3600
             wgt = (10800. - deltasec) / 10800.
          else if (UTChr(i,j) == secondUTChr) then
             TS_A = ARFS_TS_02(i,j)
             TS_B = ARFS_TS_03(i,j)
             wgt = 1
          else
             TS_A = ARFS_TS_02(i,j)
             TS_B = ARFS_TS_03(i,j)
             deltasec = ( UTChr(i,j) - secondUTChr ) * 3600
             wgt = (10800. - deltasec) / 10800.
          end if
          if (TS_A > 0 .and. TS_B > 0) then
             TS = ((wgt)*TS_A) + ((1. - wgt)*TS_B)
          else
             cycle
          end if

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

          END IF
       END DO !ii=1,nrow !ROW LON
    END DO !jj=1,mcol !COL LAT

    !write soil moisture retrieval outputs
    L1B_dir_len = len_trim(SMAPeOPL%L1Bdir)
    L1B_fname_len = len_trim(SMAPFILE)

    if(SMAPeOPL%L1Btype.eq.1) then  !NRT
       retrieval_fname = trim(SMAPeOPL%SMoutdir)//"/"//"ARFS_SM_V_"//&
                         trim(SMAPFILE(L1B_dir_len+18:L1B_fname_len-3))//".nc"
       yyyymmdd = trim(SMAPFILE(L1B_fname_len-28:L1B_fname_len-20))
       hhmmss = trim(SMAPFILE(L1B_fname_len-19:L1B_fname_len-13))
    elseif(SMAPeOPL%L1Btype.eq.2) then  !Historical
       retrieval_fname = trim(SMAPeOPL%SMoutdir)//"/"//"ARFS_SM_V_"//&
            trim(SMAPFILE(L1B_dir_len+14:L1B_fname_len-3))//".nc"
       yyyymmdd = trim(SMAPFILE(L1B_fname_len-28:L1B_fname_len-20))
       hhmmss = trim(SMAPFILE(L1B_fname_len-19:L1B_fname_len-13))
    endif

    write (LDT_logunit,*) '[INFO] Writing soil moisture retrieval file ', trim(retrieval_fname)

    ! NOTE: nrow is actually number of columns, mcol is actually number of
    ! rows
    call LDT_ARFSSM_write_netcdf(nrow, mcol, arfs_sm, retrieval_fname, &
         yyyymmdd, hhmmss)
    write (LDT_logunit,*) '[INFO] Successfully wrote soil moisture retrieval file ', trim(retrieval_fname)
    write (LDT_logunit,*) '[INFO] Finished generating soil moisture retrievals'

 end subroutine ARFSSMRETRIEVAL
