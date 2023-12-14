!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! SUBROUTINE: SMAPL1BRESAMPLE
!
! REVISION HISTORY:
!  22 Oct 2021: P.W.LIU; Initial implemetation
!  18 Dec 2021: Yonghwan Kwon; modified for LDT
!  09 Feb 2023: Eric Kemp; now processes subset of fields, no output of
!    data to separate binary files.
!
! DESCRIPTION: RESAMPLE SMAPL1B TB TO AIR FORCE GRID
! INPUT : SMAP - L1B Brightness Temperature
! OUTPUT: SMAPTB_ARFSGRIDE_ddmmyyy.dat
! NOTES : Inverse Distance Squared with 0.4 deg serching window
!-------------------------------------------------------------------------

subroutine SMAPL1BRESAMPLE_SUBSET(SMAPFILE,L1B_dir,Orbit,ARFS_TIME,rc)

  USE VARIABLES
  USE DATADOMAIN
  USE FUNCTIONS
  USE TOOLSUBS
  USE invdist_l1b2arfs
  USE LDT_logMod
  USE LDT_smap_e_oplMod

  IMPLICIT NONE

  INTEGER :: i, j, nrow, mcol
  CHARACTER (len=100) :: SMAPFILE
  character (len=100) :: L1B_dir
  character (len=20)  :: variable_name(13)
  character (len=100) :: resample_filename(13)
  character (len=1)   :: Orbit
  integer             :: var_i
  integer             :: L1B_dir_len,L1B_fname_len
  integer :: ierr
  integer :: rc

  REAL*4,DIMENSION(:,:),ALLOCATABLE :: TIME_L1B, TBV_COR_L1B
  REAL*4,DIMENSION(:,:),ALLOCATABLE :: LAT_L1B, LON_L1B, SCNANG_L1B
  REAL*4,DIMENSION(:),ALLOCATABLE :: ANTSCN_L1B
  INTEGER*4,DIMENSION(:,:),ALLOCATABLE :: TBVFLAG_L1B, TBHFLAG_L1B
  REAL*8,DIMENSION(:), ALLOCATABLE :: ARFS_LAT, ARFS_LON
  REAL*4,DIMENSION(2560,1920) :: ARFS_TIME, ARFS_COR_TBV

  REAL :: T1, T2

  rc = 0

  CALL ARFS_GEO
  ALLOCATE(ARFS_LAT(arfs_nrow_lat),ARFS_LON(arfs_mcol_lon))
  ARFS_LAT = LAT(arfs_geo_lat_lo,arfs_geo_lat_up,-arfs_lat_space)
  ARFS_LON = LON(arfs_geo_lon_lf,arfs_geo_lon_rt,arfs_lon_space)

  !Input (Path/filename, datatypes, lat, lon) Return(DATA,LAT,LON,length of row and col)
  ! READ SMAP_L1B DATA FROM HDF5

  ! EMK...Try fault tolerant NRT version.
  CALL GetSMAP_L1B_NRT_SUBSET(SMAPFILE, TIME_L1B, TBV_COR_L1B, &
       LAT_L1B, LON_L1B, TBVFLAG_L1B, TBHFLAG_L1B, &
       ANTSCN_L1B, SCNANG_L1B, nrow, mcol, ierr)
  if (ierr == 1) then
     if (nrow == 0 .and. mcol == 0) then
        write(LDT_logunit,*)'[ERR] Problem reading ', trim(SMAPFILE)
        rc = 1
        return
     else
        write(LDT_logunit,*)'[ERR] Unknown internal error!'
        write(LDT_logunit,*)'[ERR] Aborting...'
        call LDT_endrun()
     end if
  end if

  !Input (DATA,LAT,LON,length of row and col); Return(TB in ARFS GRID)
  !CALL L1BTB2ARFS_INVDIS(TIME_L1B, TBV_COR_L1B, TBH_COR_L1B, TBV_L1B, TBH_L1B, SURWAT_V_L1B, SURWAT_H_L1B, &
  !                       NETD_V_L1B, NETD_H_L1B, LAT_L1B, LON_L1B, TBVFLAG_L1B, TBHFLAG_L1B, ANTSCN_L1B, SCNANG_L1B, nrow, mcol, &
  !                       ARFS_LAT, ARFS_LON, ARFS_TIME, ARFS_COR_TBV, ARFS_COR_TBH, ARFS_TBV, ARFS_TBH, ARFS_NEDTV, ARFS_NEDTH, &
  !                       ARFS_SURWAT_V, ARFS_SURWAT_H, ARFS_WTV, ARFS_WTH, ARFS_SAMPLE_V, ARFS_SAMPLE_H)
  CALL L1BTB2ARFS_INVDIS_SUBSET(TIME_L1B, TBV_COR_L1B, &
       LAT_L1B, LON_L1B, TBVFLAG_L1B, TBHFLAG_L1B, ANTSCN_L1B, SCNANG_L1B, &
       nrow, mcol, &
       ARFS_LAT, ARFS_LON, ARFS_TIME, ARFS_COR_TBV)
  SMAPeOPL%ARFS_TBV_COR = ARFS_COR_TBV

  L1B_dir_len = len_trim(L1B_dir)
  L1B_fname_len = len_trim(SMAPFILE)

  if(SMAPeOPL%L1Btype.eq.1) then  !NRT
     Orbit = trim(SMAPFILE(L1B_dir_len+24:L1B_dir_len+24))
  elseif(SMAPeOPL%L1Btype.eq.2) then  !Historical
     Orbit = trim(SMAPFILE(L1B_dir_len+20:L1B_dir_len+20))
  endif

  ! Cleanup
  deallocate(TBV_COR_L1B)
  deallocate(LAT_L1B)
  deallocate(LON_L1B)
  deallocate(SCNANG_L1B)
  deallocate(ANTSCN_L1B)
  deallocate(TBVFLAG_L1B)
  deallocate(TBHFLAG_L1B)
  deallocate(ARFS_LAT)
  deallocate(ARFS_LON)

end subroutine SMAPL1BRESAMPLE_SUBSET
