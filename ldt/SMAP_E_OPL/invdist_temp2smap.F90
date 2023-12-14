!=======================================================================
!  MODULE, invdist_temp2smap P.W.LIU, 01/11/2022
!  SUBROUTINE TO RESAMPLE SOIL TEMPERATURE FROM NOAH 39 ON ARFS GRID TO
!  SMAP RESOLUTIN AT ARFS GRID
!=======================================================================
 MODULE invdist_temp2smap
   IMPLICIT NONE

 CONTAINS
   !SUBROUTINE RESAMPLETEMP(arfs_temp,arfs_lat,arfs_lon,arfs_fine_lat,arfs_fine_lon,arfs_fine_wt,arfs_fine_temp,arfs_wt,arfs_smap_temp)
   SUBROUTINE RESAMPLETEMP(arfs_temp,arfs_lat,arfs_lon,arfs_fine_lat,arfs_fine_lon,arfs_smap_temp)
 

     INTEGER(4) :: ii, jj, k, r, c, rr, rmin, rmax, cc, cmin, cmax, nrows_l1btb, ncols_l1btb
     REAL(8), PARAMETER :: RE_KM = 6371.228, search_radius = 20.0, PI = 3.141592653589793238, d2r = PI/180.0
     REAL(8)  :: gcdist, lat1, lon1, lat2, lon2
     REAL*8 ,DIMENSION(:), ALLOCATABLE :: arfs_lat, arfs_lon
     REAL*8 ,DIMENSION(:), ALLOCATABLE :: arfs_fine_lat, arfs_fine_lon
     REAL*4,DIMENSION(2560,1920) :: arfs_temp, arfs_smap_temp, arfs_wt
     REAL*4,DIMENSION(7680,5760) :: arfs_fine_temp, arfs_fine_wt
     INTEGER(4),DIMENSION(2560,1920):: zerodistflag
     INTEGER(4),DIMENSION(7680,5760):: zerodistflag_fine

     !INITIAL THE OUTPUT VARIABLES
     zerodistflag=0
     zerodistflag_fine=0
     arfs_smap_temp=0
     arfs_wt=0
     arfs_fine_temp=0
     arfs_fine_wt=0
     !DISAGGREGATE ARFS TEMP TO FINEER GRID (~3km)
     DO ii=1,1920
        !DO ii=522,524
        lat1 = DBLE (arfs_lat(ii)*d2r)
        c = MINLOC(ABS(arfs_lat(ii)-arfs_fine_lat(:)),1) !Lat Direction r: row in fine scale
        !PRINT*, arfs_fine_lat(r), arfs_lat(ii)
        cmin=c-3 ; IF (cmin < 1) cmin=1
        cmax=c+3 ; IF (cmax > size(arfs_fine_lat)) cmax=size(arfs_fine_lat)

        DO jj=1,2560
           !DO jj=82,84
           lon1 = DBLE (arfs_lon(jj)*d2r)
           r = MINLOC(ABS(arfs_lon(jj)-arfs_fine_lon(:)),1) !Lat Direction c: colume in fine scale
           !PRINT*, arfs_fine_lon(c), arfs_lon(jj)
           rmin=r-3 ; IF (rmin < 1) rmin=1
           rmax=r+3 ; IF (rmax > size(arfs_fine_lon)) rmax=size(arfs_fine_lon)
           !PRINT*,'cmin cmax', cmin, cmax, size(arfs_fine_lon)
           !DO rr=rmin,rmax
           !   DO cc=cmin,cmax
           DO cc=cmin,cmax
              lat2 = DBLE (arfs_fine_lat(cc)*d2r)
              DO rr=rmin,rmax
                 lon2 = DBLE (arfs_fine_lon(rr)*d2r)
                 !PRINT*, lat1, lon1, lat2, lon2
                 !gcdist = RE_KM * DACOS ( DSIN (lat1) * DSIN (lat2) + DCOS (lat1) * DCOS (lat2) * DCOS (lon1-lon2) )  !original
                 !---------------------------------------kyh
                 if(lat1.eq.lat2.and.lon1.eq.lon2) then
                    gcdist = 0.
                 else
                    gcdist = RE_KM * DACOS ( DSIN (lat1) * DSIN (lat2) + DCOS (lat1) * DCOS (lat2) * DCOS (lon1-lon2) )
                 endif
                 !---------------------------------------kyh
                 !PRINT*, lat1, lon1, lat2, lon2, gcdist
                 !PRINT*, arfs_temp(ii,jj)
                 IF (gcdist < search_radius) THEN !RESAMPLE ONLY WITHIN THE SEARCH RANGE
                    !PRINT*,'HERE 0'
                    IF (gcdist < 0.0001D0) THEN !The TB is right on the grid center
                       zerodistflag_fine(rr,cc) = 1
                       !PRINT*,'HERE 01'
                       IF ((ABS (arfs_temp(jj,ii)-(-9999.000)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(0)
                          arfs_fine_temp(rr,cc) = arfs_temp(jj,ii) ; arfs_fine_wt(rr,cc) = 1.0
                          !PRINT*,'Here 1'
                       ENDIF
                    ELSE
                       IF (zerodistflag_fine(rr,cc).EQ.0) THEN !To maintain the corresponding pixel has the same value
                          IF ((ABS (arfs_temp(jj,ii)-(-9999.000)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(0)
                             arfs_fine_temp(rr,cc) = arfs_fine_temp(rr,cc)+arfs_temp(jj,ii) / SNGL(gcdist)
                             arfs_fine_wt(rr,cc) = arfs_fine_wt(rr,cc) + 1.0 / SNGL(gcdist)
                             !PRINT*,'Here 2'
                          ENDIF
                       ENDIF !(zerodistflag_fine(rr,cc).NE.0)
                    ENDIF !IF (gcdist < 0.0001D0)
                 ENDIF !(gcdist < search_radius)
              ENDDO !rr=rmin,rmax
           ENDDO !cc=cmin,cmax
        ENDDO !jj=1:2560
     ENDDO !ii=1,1920
     WHERE(arfs_fine_temp.NE.0.0.AND.arfs_fine_wt.NE.0.0)
        arfs_fine_temp=arfs_fine_temp/arfs_fine_wt
     ENDWHERE
     !UPSCALL ARFS FINE TEMP TO ARFS GRID (at ~33km)
     DO ii=1,1920
        c=3*ii-1  !3x3 FINE GRID EQ TO A ARFS GRID
        cmin=c-5 ; IF (cmin < 1) cmin=1
        cmax=c+5 ; IF (cmax > size(arfs_fine_lat)) cmax=size(arfs_fine_lat)
        DO jj=1,2560
           r=3*jj-1  !3x3 FINE GRID EQ TO A ARFS GRID
           rmin=r-5 ; IF (rmin < 1) rmin=1
           rmax=r+5 ; IF (rmax > size(arfs_fine_lon)) rmax=size(arfs_fine_lon)
           DO cc=cmin,cmax
              DO rr=rmin,rmax
                 IF (abs(arfs_fine_temp(rr,cc)).GT.1.0D-7) THEN !DO WHEN T ~0 (change to fillvalue)
                    arfs_smap_temp(jj,ii)=arfs_smap_temp(jj,ii)+arfs_fine_temp(rr,cc); arfs_wt(jj,ii)=arfs_wt(jj,ii)+1.0 !Weight the point 1
                 ENDIF
              ENDDO !rr=rmin,rmax
           ENDDO !cc=cmin,cmax
        ENDDO !jj=1:2560
     ENDDO !ii=1:1920
     WHERE(arfs_smap_temp.NE.0.0.AND.arfs_wt.NE.0.0)
        arfs_smap_temp=arfs_smap_temp/arfs_wt
     ENDWHERE
     WHERE(arfs_smap_temp.EQ.0.0)
        arfs_smap_temp=-9999
     ENDWHERE

   END SUBROUTINE RESAMPLETEMP
 END MODULE invdist_temp2smap
