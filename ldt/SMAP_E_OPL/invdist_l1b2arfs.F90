!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT)
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: invdist_l1b2arfs
!
! REVISION HISTORY: 
!  26 Oct 2021: P.W.LIU; Initial implemetation
!  28 Jan 2022: P.W.LIU; CHNAGE OUTPUT COORDINATE TO COMPLY LIS OUTPUT
!  22 Feb 2022: Yonghwan Kwon; modified for LDT
!  10 Feb 2023: Eric Kemp; Added code to resample subset of variables.
!
! DESCRIPTION: SUBROUTINE TO RESAMPLE L1B TB ONTO AIR FORCE GRID
!-------------------------------------------------------------------------

 MODULE invdist_l1b2arfs
   IMPLICIT NONE

 CONTAINS
   SUBROUTINE L1BTB2ARFS_INVDIS(tim, tbvl1b_cor, tbhl1b_cor, tbvl1b, tbhl1b, surwat_v_l1b, surwat_h_l1b, &
        netd_v_l1b, netd_h_l1b, lat_l1b, lon_l1b, tbv_qual_flag, tbh_qual_flag, sc_nadir_angle, antenna_scan_angle, nrows_l1btb, ncols_l1btb, &
        ref_lat, ref_lon, arfs_tim, arfs_tbv_cor, arfs_tbh_cor, arfs_tbv, arfs_tbh, arfs_nedtv, arfs_nedth, &
        arfs_surwatv, arfs_surwath, arfs_wt_cor_tbv, arfs_wt_cor_tbh, arfs_samplenumv, arfs_samplenumh)


     INTEGER(4) :: ii, jj, k, r, c, rr, rmin, rmax, cc, cmin, cmax, nrows_l1btb, ncols_l1btb
     INTEGER(4), PARAMETER :: qualitybit = 0
     REAL(8), PARAMETER :: RE_KM = 6371.228, search_radius = 20.0, PI = 3.141592653589793238, d2r = PI/180.0
     REAL(8)  :: gcdist, lat1, lon1, lat2, lon2
     REAL*4,DIMENSION(nrows_l1btb,ncols_l1btb) :: tim, tbvl1b_cor, tbhl1b_cor, tbvl1b, tbhl1b, surwat_v_l1b, surwat_h_l1b, netd_v_l1b, netd_h_l1b
     REAL*4,DIMENSION(nrows_l1btb,ncols_l1btb) :: lat_l1b, lon_l1b, antenna_scan_angle
     REAL*4,DIMENSION(ncols_l1btb) :: sc_nadir_angle
     INTEGER*4,DIMENSION(nrows_l1btb,ncols_l1btb) :: tbv_qual_flag, tbh_qual_flag
     INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: zerodistflag
     REAL*8,DIMENSION(:),ALLOCATABLE :: ref_lat, ref_lon
     REAL*4,DIMENSION(2560,1920) :: arfs_tim, arfs_tbv_cor, arfs_tbh_cor, arfs_tbv, arfs_tbh, arfs_nedtv, arfs_nedth, arfs_surwatv, arfs_surwath
     REAL*4,DIMENSION(2560,1920) :: arfs_wt_tim, arfs_wt_cor_tbv, arfs_wt_cor_tbh, arfs_wt_tbv, arfs_wt_tbh, arfs_wt_nedtv, arfs_wt_nedth, arfs_wt_surwatv, arfs_wt_surwath
     INTEGER*4,DIMENSION(2560,1920) :: arfs_samplenumv, arfs_samplenumh

     !ALLOCATE(zerodistflag(size(ref_lat),size(ref_lon)))
     ALLOCATE(zerodistflag(size(ref_lon),size(ref_lat)))
     !INITIAL THE OUTPUT VARIABLES
     arfs_tim=0.0
     arfs_tbv_cor=0.0
     arfs_tbh_cor=0.0
     arfs_tbv=0.0
     arfs_tbh=0.0
     arfs_nedtv=0.0
     arfs_nedth=0.0
     arfs_surwatv=0.0
     arfs_surwath=0.0
     arfs_wt_tim=0.0
     arfs_wt_cor_tbv=0.0
     arfs_wt_cor_tbh=0.0
     arfs_wt_tbv=0.0
     arfs_wt_tbh=0.0
     arfs_wt_nedtv=0.0
     arfs_wt_nedth=0.0
     arfs_wt_surwatv=0.0
     arfs_wt_surwath=0.0
     arfs_samplenumv=0.0

     DO ii = 1,ncols_l1btb
        IF (ABS (sc_nadir_angle(ii)) <= 2.0) THEN
           DO jj = 1,nrows_l1btb
              IF (ABS (antenna_scan_angle(jj,ii)).LE.360.00) THEN
                 ! FIND ARFS_GRID (r,c)
                 c = MINLOC(ABS(lat_l1b(jj,ii)-ref_lat(:)),1) !Lat Direction
                 r = MINLOC(ABS(lon_l1b(jj,ii)-ref_lon(:)),1) !Lon Direction
                 rmin=r-5 ; IF (rmin < 1) rmin=1
                 rmax=r+5 ; IF (rmax > size(ref_lon)) rmax=size(ref_lon)
                 cmin=c-5 ; IF (cmin < 1) cmin=1
                 cmax=c+5 ; IF (cmax > size(ref_lat)) cmax=size(ref_lat)
                 IF (IBITS (tbv_qual_flag(jj,ii),qualitybit,1) == 0 .AND. IBITS (tbh_qual_flag(jj,jj),qualitybit,1) == 0) THEN !RESAMPLE ONLY WHEN BOTH V and H MEET QUALITY
                    k=0
                    DO rr = rmin,rmax !Lon direction
                       DO cc =cmin,cmax !Lat direction
                          lat1 = DBLE (lat_l1b(jj,ii)*d2r)
                          lon1 = DBLE (lon_l1b(jj,ii)*d2r)
                          lat2 = DBLE (ref_lat(cc)*d2r)
                          lon2 = DBLE (ref_lon(rr)*d2r)

                          if(lat1.eq.lat2.and.lon1.eq.lon2) then
                             gcdist = 0.
                          else
                             gcdist = RE_KM * DACOS ( DSIN (lat1) * DSIN (lat2) + DCOS (lat1) * DCOS (lat2) * DCOS (lon1-lon2) )
                          endif

                          IF (gcdist < search_radius) THEN !RESAMPLE ONLY WITHIN THE SEARCH RANGE
                             IF (gcdist < 0.0001D0) THEN !The TB is right on the grid center
                                zerodistflag (rr,cc) = 1
                                IF ((ABS (tim(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_tim(rr,cc) = tim(jj,ii) ; arfs_wt_tim(rr,cc) = 1.0
                                END IF
                                IF ((ABS (surwat_v_l1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_surwatv(rr,cc) = surwat_v_l1b(jj,ii) ; arfs_wt_surwatv(rr,cc) = 1.0
                                END IF
                                IF ((ABS (surwat_h_l1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_surwath(rr,cc) = surwat_h_l1b(jj,ii) ; arfs_wt_surwath(rr,cc) = 1.0
                                END IF
                                IF ((ABS (netd_v_l1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_nedtv(rr,cc) = netd_v_l1b(jj,ii) ; arfs_wt_nedtv(rr,cc) = 1.0
                                END IF
                                IF ((ABS (netd_h_l1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_nedth(rr,cc) = netd_h_l1b(jj,ii) ; arfs_wt_nedth(rr,cc) = 1.0
                                END IF
                                IF ((ABS (tbvl1b_cor(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_tbv_cor(rr,cc) = tbvl1b_cor(jj,ii) ; arfs_wt_cor_tbv(rr,cc) = 1.0
                                   arfs_samplenumv(rr,cc)=1 !Sample number only calculate for correct tb
                                   k=k+1;
                                END IF
                                IF ((ABS (tbhl1b_cor(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_tbh_cor(rr,cc) = tbhl1b_cor(jj,ii) ; arfs_wt_cor_tbh(rr,cc) = 1.0
                                   arfs_samplenumh(rr,cc)=1 !Sample number only calculate for correct tb
                                END IF
                                IF ((ABS (tbvl1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_tbv(rr,cc) = tbvl1b(jj,ii) ; arfs_wt_tbv(rr,cc) = 1.0
                                END IF
                                IF ((ABS (tbhl1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_tbh(rr,cc) = tbhl1b(jj,ii) ; arfs_wt_tbh(rr,cc) = 1.0
                                END IF
                             ELSE
                                IF (zerodistflag (rr,cc).EQ.0) THEN

                                   IF ((ABS (tim(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_tim(rr,cc) = arfs_tim(rr,cc) + tim(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_tim(rr,cc) = arfs_wt_tim(rr,cc) + 1.0 / SNGL (gcdist*gcdist)
                                   END IF
                                   IF ((ABS (surwat_v_l1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_surwatv(rr,cc) = arfs_surwatv(rr,cc) + surwat_v_l1b(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_surwatv(rr,cc) = arfs_wt_surwatv(rr,cc) + 1.0 / SNGL (gcdist*gcdist)
                                   END IF
                                   IF ((ABS (surwat_h_l1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_surwath(rr,cc) = arfs_surwath(rr,cc) + surwat_h_l1b(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_surwath(rr,cc) = arfs_wt_surwath(rr,cc) + 1.0 / SNGL (gcdist*gcdist)
                                   END IF
                                   IF ((ABS (netd_v_l1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_nedtv(rr,cc) = arfs_nedtv(rr,cc) + netd_v_l1b(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_nedtv(rr,cc) = arfs_wt_nedtv(rr,cc) + 1.0 / SNGL (gcdist*gcdist)
                                   END IF
                                   IF ((ABS (netd_h_l1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_nedth(rr,cc) = arfs_nedth(rr,cc) + netd_h_l1b(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_nedth(rr,cc) = arfs_wt_nedth(rr,cc) + 1.0 / SNGL (gcdist*gcdist)
                                   END IF
                                   IF ((ABS (tbvl1b_cor(jj,ii) - (-9999.0))).GT.1.0D-7) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_tbv_cor(rr,cc) = arfs_tbv_cor(rr,cc) + tbvl1b_cor(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_cor_tbv(rr,cc) = arfs_wt_cor_tbv(rr,cc) + 1.0 /  SNGL (gcdist*gcdist)
                                      arfs_samplenumv(rr,cc)=arfs_samplenumv(rr,cc)+1.0 !Sample number only calculate for correct tb
                                      k=k+1;
                                   END IF
                                   IF ((ABS (tbhl1b_cor(jj,ii) - (-9999.0))).GT.1.0D-7) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_tbh_cor(rr,cc) = arfs_tbh_cor(rr,cc) + tbhl1b_cor(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_cor_tbh(rr,cc) = arfs_wt_cor_tbh(rr,cc) + 1.0 /  SNGL (gcdist*gcdist)
                                      arfs_samplenumh(rr,cc)=arfs_samplenumh(rr,cc)+1.0 !Sample number only calculate for correct tb
                                   END IF
                                   IF ((ABS (tbvl1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_tbv(rr,cc) = arfs_tbv(rr,cc) + tbvl1b(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_tbv(rr,cc) = arfs_wt_tbv(rr,cc) + 1.0 / SNGL (gcdist*gcdist)
                                   END IF
                                   IF ((ABS (tbhl1b(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_tbh(rr,cc) = arfs_tbh(rr,cc) + tbhl1b(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_tbh(rr,cc) = arfs_wt_tbh(rr,cc) + 1.0 / SNGL (gcdist*gcdist)
                                   END IF
                                END IF !(zerodistflag (rr,cc) = 0)
                             END IF !(gcdist < 0.0001D0)!
                          END IF !(gcdist < search_radius)
                       END DO !cc =cmin,cmax
                    END DO !rr = rmin,rmax
                 END IF !(IBITS (tbv_qual_flag(jj,ii),qualitybit,1) == 0 .AND. IBITS (tbh_qual_flag(ii,jj),qualitybit,1) == 0)
              END IF !(ABS (antenna_scan_angle(ii,jj)) <= 360.00)
           END DO !jj=1,2
        END IF !(ABS (sc_nadir_angle(ii)) <= 2.0)
     END DO !ii=1,2

     !APPLY WEIGHTING FUNCTION FOR RESAMPLING
     WHERE(arfs_tim.NE.0.0.AND.arfs_wt_tim.NE.0.0)
        arfs_tim = arfs_tim / arfs_wt_tim
     END WHERE
     WHERE(arfs_surwatv.NE.0.0.AND.arfs_wt_surwatv.NE.0.0)
        arfs_surwatv = arfs_surwatv / arfs_wt_surwatv
     END WHERE
     WHERE(arfs_surwath.NE.0.0.AND.arfs_wt_surwath.NE.0.0)
        arfs_surwath = arfs_surwath / arfs_wt_surwath
     END WHERE
     WHERE(arfs_nedtv.NE.0.0.AND.arfs_wt_nedtv.NE.0.0)
        arfs_nedtv = arfs_nedtv / arfs_wt_nedtv
     END WHERE
     WHERE(arfs_nedth.NE.0.0.AND.arfs_wt_nedth.NE.0.0)
        arfs_nedth = arfs_nedth / arfs_wt_nedth
     END WHERE
     WHERE(arfs_tbv_cor.NE.0.0.AND.arfs_wt_cor_tbv.NE.0.0)
        arfs_tbv_cor = arfs_tbv_cor / arfs_wt_cor_tbv
     END WHERE
     WHERE(arfs_tbh_cor.NE.0.0.AND.arfs_wt_cor_tbh.NE.0.0)
        arfs_tbh_cor = arfs_tbh_cor / arfs_wt_cor_tbh
     END WHERE
     WHERE(arfs_tbv.NE.0.0.AND.arfs_wt_tbv.NE.0.0)
        arfs_tbv = arfs_tbv / arfs_wt_tbv
     END WHERE
     WHERE(arfs_tbh.NE.0.0.AND.arfs_wt_tbh.NE.0.0)
        arfs_tbh = arfs_tbh / arfs_wt_tbh
     END WHERE

   END SUBROUTINE L1BTB2ARFS_INVDIS

   ! EMK...Only process subset of SMAP L1B fields for NRT operations.
   SUBROUTINE L1BTB2ARFS_INVDIS_SUBSET(tim, tbvl1b_cor, &
        lat_l1b, lon_l1b, tbv_qual_flag, tbh_qual_flag, &
        sc_nadir_angle, antenna_scan_angle, nrows_l1btb, ncols_l1btb, &
        ref_lat, ref_lon, arfs_tim, arfs_tbv_cor)

     INTEGER(4) :: ii, jj, k, r, c, rr, rmin, rmax, cc, cmin, cmax, nrows_l1btb, ncols_l1btb
     INTEGER(4), PARAMETER :: qualitybit = 0
     REAL(8), PARAMETER :: RE_KM = 6371.228, search_radius = 20.0, PI = 3.141592653589793238, d2r = PI/180.0
     REAL(8)  :: gcdist, lat1, lon1, lat2, lon2
     REAL*4,DIMENSION(nrows_l1btb,ncols_l1btb) :: tim, tbvl1b_cor
     REAL*4,DIMENSION(nrows_l1btb,ncols_l1btb) :: lat_l1b, lon_l1b, antenna_scan_angle
     REAL*4,DIMENSION(ncols_l1btb) :: sc_nadir_angle
     INTEGER*4,DIMENSION(nrows_l1btb,ncols_l1btb) :: tbv_qual_flag, tbh_qual_flag
     INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: zerodistflag
     REAL*8,DIMENSION(:),ALLOCATABLE :: ref_lat, ref_lon
     REAL*4,DIMENSION(2560,1920) :: arfs_tim, arfs_tbv_cor
     REAL*4,DIMENSION(2560,1920) :: arfs_wt_tim, arfs_wt_cor_tbv
     INTEGER*4,DIMENSION(2560,1920) :: arfs_samplenumv

     ALLOCATE(zerodistflag(size(ref_lon),size(ref_lat)))
     !INITIAL THE OUTPUT VARIABLES
     arfs_tim=0.0
     arfs_tbv_cor=0.0
     arfs_wt_tim=0.0
     arfs_wt_cor_tbv=0.0
     arfs_samplenumv=0.0

     DO ii = 1,ncols_l1btb
        IF (ABS (sc_nadir_angle(ii)) <= 2.0) THEN
           DO jj = 1,nrows_l1btb
              IF (ABS (antenna_scan_angle(jj,ii)).LE.360.00) THEN
                 lat1 = DBLE (lat_l1b(jj,ii)*d2r)
                 lon1 = DBLE (lon_l1b(jj,ii)*d2r)
                 ! FIND ARFS_GRID (r,c)
                 c = MINLOC(ABS(lat_l1b(jj,ii)-ref_lat(:)),1) !Lat Direction
                 r = MINLOC(ABS(lon_l1b(jj,ii)-ref_lon(:)),1) !Lon Direction
                 rmin=r-5 ; IF (rmin < 1) rmin=1
                 rmax=r+5 ; IF (rmax > size(ref_lon)) rmax=size(ref_lon)
                 cmin=c-5 ; IF (cmin < 1) cmin=1
                 cmax=c+5 ; IF (cmax > size(ref_lat)) cmax=size(ref_lat)
!                 IF (IBITS (tbv_qual_flag(jj,ii),qualitybit,1) == 0 .AND. IBITS (tbh_qual_flag(jj,jj),qualitybit,1) == 0) THEN !RESAMPLE ONLY WHEN BOTH V and H MEET QUALITY
                 IF (IBITS (tbv_qual_flag(jj,ii),qualitybit,1) == 0 .AND. IBITS (tbh_qual_flag(jj,ii),qualitybit,1) == 0) THEN !RESAMPLE ONLY WHEN BOTH V and H MEET QUALITY
                    k=0
                    !                         DO rr = rmin,rmax !Lon direction
                    !                            DO cc =cmin,cmax !Lat direction
                    DO cc =cmin,cmax !Lat direction
                       lat2 = DBLE (ref_lat(cc)*d2r)
                       DO rr = rmin,rmax !Lon direction
                          lon2 = DBLE (ref_lon(rr)*d2r)

                          if(lat1.eq.lat2.and.lon1.eq.lon2) then
                             gcdist = 0.
                          else
                             gcdist = RE_KM * DACOS ( DSIN (lat1) * DSIN (lat2) + DCOS (lat1) * DCOS (lat2) * DCOS (lon1-lon2) )
                          endif

                          IF (gcdist < search_radius) THEN !RESAMPLE ONLY WITHIN THE SEARCH RANGE
                             IF (gcdist < 0.0001D0) THEN !The TB is right on the grid center
                                zerodistflag (rr,cc) = 1
                                IF ((ABS (tim(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_tim(rr,cc) = tim(jj,ii) ; arfs_wt_tim(rr,cc) = 1.0
                                END IF
                                IF ((ABS (tbvl1b_cor(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   arfs_tbv_cor(rr,cc) = tbvl1b_cor(jj,ii) ; arfs_wt_cor_tbv(rr,cc) = 1.0
                                   arfs_samplenumv(rr,cc)=1 !Sample number only calculate for correct tb
                                   k=k+1;
                                END IF
                             ELSE
                                IF (zerodistflag (rr,cc).EQ.0) THEN

                                   IF ((ABS (tim(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_tim(rr,cc) = arfs_tim(rr,cc) + tim(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_tim(rr,cc) = arfs_wt_tim(rr,cc) + 1.0 / SNGL (gcdist*gcdist)
                                   END IF
                                   IF ((ABS (tbvl1b_cor(jj,ii) - (-9999.0))).GT.1.0D-7) THEN !DO IF NOT FILLVALUE(-9999)
                                      arfs_tbv_cor(rr,cc) = arfs_tbv_cor(rr,cc) + tbvl1b_cor(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_cor_tbv(rr,cc) = arfs_wt_cor_tbv(rr,cc) + 1.0 /  SNGL (gcdist*gcdist)
                                      arfs_samplenumv(rr,cc)=arfs_samplenumv(rr,cc)+1.0 !Sample number only calculate for correct tb
                                      k=k+1;
                                   END IF

                                END IF !(zerodistflag (rr,cc) = 0)
                             END IF !(gcdist < 0.0001D0)!
                          END IF !(gcdist < search_radius)
                          !                            END DO !cc =cmin,cmax
                          !                         END DO !rr = rmin,rmax
                       END DO !rr = rmin,rmax
                    END DO !cc =cmin,cmax

                 END IF !(IBITS (tbv_qual_flag(jj,ii),qualitybit,1) == 0 .AND. IBITS (tbh_qual_flag(ii,jj),qualitybit,1) == 0)
              END IF !(ABS (antenna_scan_angle(ii,jj)) <= 360.00)
           END DO !jj=1,2
        END IF !(ABS (sc_nadir_angle(ii)) <= 2.0)
     END DO !ii=1,2

     !APPLY WEIGHTING FUNCTION FOR RESAMPLING
     WHERE(arfs_tim.NE.0.0.AND.arfs_wt_tim.NE.0.0)
        arfs_tim = arfs_tim / arfs_wt_tim
     END WHERE
     WHERE(arfs_tbv_cor.NE.0.0.AND.arfs_wt_cor_tbv.NE.0.0)
        arfs_tbv_cor = arfs_tbv_cor / arfs_wt_cor_tbv
     END WHERE

     deallocate(zerodistflag) ! EMK cleanup memory
   END SUBROUTINE L1BTB2ARFS_INVDIS_SUBSET

 END MODULE invdist_l1b2arfs
