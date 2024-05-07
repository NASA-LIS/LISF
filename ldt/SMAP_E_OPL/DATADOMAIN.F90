!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
MODULE DATADOMAIN
       USE VARIABLES
       IMPLICIT NONE
CONTAINS

     SUBROUTINE GLDAS025d_DOMAIN
     !DOMAIN OF GLDAS 0.25d DATA
     ldas_lon_lf = -180
     ldas_lon_rt = 180
     ldas_lat_up = 90
     ldas_lat_lo = -60
     ldas_space = 0.25
     row_ldas= 1440
     col_ldas= 600
     RETURN
     !------------------------------
     END SUBROUTINE GLDAS025d_DOMAIN


     SUBROUTINE NLDAS_DOMAIN
     !DOMAIN OF NLDAS DATA
     ldas_lon_lf = -125
     ldas_lon_rt = -67
     ldas_lat_up = 53
     ldas_lat_lo = 25
     ldas_space = 0.125
     row_ldas= 464
     col_ldas= 224
     RETURN
     !------------------------------
     END SUBROUTINE NLDAS_DOMAIN

     SUBROUTINE AVHR_DOMAIN
     !DOMAIN OF AVHRR-NDVI DATA
      avhr_lon_lf = -180
      avhr_lon_rt = 180
      avhr_lat_up = 90
      avhr_lat_lo = -90
      avhr_space=0.05
      row_avhr = 7200
      col_avhr = 3600
     !------------------------------
     END SUBROUTINE AVHR_DOMAIN

     SUBROUTINE MODIS1KM_DOMAIN
     !DOMAIN OF MODIS1KM GLOBAL DATA
!      modis_lon_lf = -180 !Global
!      modis_lon_rt = 180
!      modis_lat_up = 90
!      modis_lat_lo = -90
!      modis_space=0.01
!      row_modis = 36000
!      col_modis = 18000

     !DOMAIN OF MODIS1KM CONUS DATA
      modis_lon_lf = -126 !NLDAS
      modis_lon_rt = -66
      modis_lat_up = 54
      modis_lat_lo = 24
      modis_space=0.01
      row_modis = 6000
      col_modis = 3000

!       modis_lon_lf = -160.5 !HAWAII
!       modis_lon_rt = -154.5
!       modis_lat_up = 23.5
!       modis_lat_lo = 18.5
!       modis_space=0.01
!       row_modis = 600
!       col_modis = 500
     !------------------------------
     END SUBROUTINE MODIS1KM_DOMAIN

     SUBROUTINE SMAP36KM_DOMAIN
     INTEGER :: i
      smap_row = 964
      smap_col = 406
      ALLOCATE(smap_lat(smap_col),smap_lon(smap_row))
      !READ Coordinates of EASE GRIDS 2,0
      OPEN(UNIT=1,FILE="DATA/EASEGRIDS/EASE2_M36km.lats.964x406x1.double", FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8,STATUS='OLD')
      DO i=1, smap_col
         READ(1,rec=(i-1)*smap_row+1) smap_lat(i)
      END DO
      CLOSE(1)

      OPEN(UNIT=1,FILE="DATA/EASEGRIDS/EASE2_M36km.lons.964x406x1.double", FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8,STATUS='OLD')
      DO i=1, smap_row
         READ(1,rec=i) smap_lon(i)
      END DO
      CLOSE(1)
      
      smap_space = 0.36
     END SUBROUTINE SMAP36KM_DOMAIN

     SUBROUTINE SMAP9KM_DOMAIN
     INTEGER :: i
      smap_row = 3856
      smap_col = 1624
      ALLOCATE(smap_lat(smap_col),smap_lon(smap_row))
      !READ Coordinates of EASE GRIDS 2,0
      OPEN(UNIT=1,FILE="DATA/EASEGRIDS/EASE2_M09km.lats.3856x1624x1.double", FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8,STATUS='OLD')
      DO i=1, smap_col
         READ(1,rec=(i-1)*smap_row+1) smap_lat(i)
      END DO
      CLOSE(1)

      OPEN(UNIT=1,FILE="DATA/EASEGRIDS/EASE2_M09km.lons.3856x1624x1.double", FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8,STATUS='OLD')
      DO i=1, smap_row
         READ(1,rec=i) smap_lon(i)
      END DO
      CLOSE(1)
      smap_space = 0.33
     END SUBROUTINE SMAP9KM_DOMAIN

     SUBROUTINE SMAP9KM_GEO
     !DOMAIN OF MODIS1KM DATA
      smap_geo_lon_lf = -180
      smap_geo_lon_rt = 180
      smap_geo_lat_up = 90
      smap_geo_lat_lo = -90
      smap_geo_space=0.10
      smap_geo_row = 3600
      smap_geo_col = 1800
     !------------------------------
     END SUBROUTINE SMAP9KM_GEO

     SUBROUTINE ARFS_GEO
     !DOMAIN OF AIRFORCE PROJECT
      arfs_geo_lon_lf = -180
      arfs_geo_lon_rt = 180
      arfs_geo_lat_up = 90
      arfs_geo_lat_lo = -90
      arfs_lon_space = 0.1406250
      arfs_lat_space = 0.0937500
      arfs_nrow_lat = 1920
      arfs_mcol_lon = 2560
     !------------------------------
     END SUBROUTINE ARFS_GEO

     SUBROUTINE ARFS_3KM_GEO
     !DOMAIN OF AIRFORCE PROJECT
      arfs_geo_lon_lf = -180
      arfs_geo_lon_rt = 180
      arfs_geo_lat_up = 90
      arfs_geo_lat_lo = -90
      arfs_lon_3km_space = 0.0468750
      arfs_lat_3km_space = 0.0312500
      arfs_nrow_lat_3km = 5760
      arfs_mcol_lon_3km = 7680
     !------------------------------
     END SUBROUTINE ARFS_3KM_GEO

     SUBROUTINE EASE1KM
     INTEGER :: i
     row_ease_1km = 34704
     col_ease_1km = 14616
      ALLOCATE(EASE1KM_lat(col_ease_1km),EASE1KM_lon(row_ease_1km))
      !READ Coordinates of EASE GRIDS 2,0
      OPEN(UNIT=1,FILE="DATA/EASEGRIDS/EASE2_M01km.lats.34704x14616x1.double", FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8,STATUS='OLD')
      DO i=1, col_ease_1km
         READ(1,rec=(i-1)*row_ease_1km+1) EASE1KM_lat(i)
      END DO
      CLOSE(1)

      OPEN(UNIT=1,FILE="DATA/EASEGRIDS/EASE2_M01km.lons.34704x14616x1.double", FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8,STATUS='OLD')
      DO i=1, row_ease_1km
         READ(1,rec=i) EASE1KM_lon(i)
      END DO
      CLOSE(1)
      END SUBROUTINE EASE1KM

     SUBROUTINE ROI_DOMAIN
     !DOMAIN OF REGION OF INTEREST
!     ROI_lon_lf = -99  !Little Washita
!     ROI_lon_rt = -97
!     ROI_lat_up = 36
!     ROI_lat_lo = 34
 
     ROI_lon_lf = -125 !NLDAS Domain
     ROI_lon_rt = -67
     ROI_lat_up = 53
     ROI_lat_lo = 25

!     ROI_lon_lf = -160.5 !Hawaii Domain
!     ROI_lon_rt = -154.5
!     ROI_lat_up = 23.5
!     ROI_lat_lo = 18.5
!     ROI_lon_lf = -112.5 !2nd_TimeZone
!     ROI_lon_rt = -97.5
!     ROI_lat_up = 53
!     ROI_lat_lo = 25

     END SUBROUTINE ROI_DOMAIN
     
     SUBROUTINE IND_1deg_1km
     INTEGER :: i     
      row_modis = 36000
      col_modis = 18000
      ALLOCATE(IND_1km_ease_lat(col_modis),IND_1km_ease_lon(row_modis))
      OPEN(UNIT=1,FILE="DATA/MODIS/IND_1deg_1km_EASE_lon.dat",STATUS='OLD')

      DO i=1, row_modis
         READ(1,*) IND_1km_ease_lon(i)
      END DO
      CLOSE(1)
 
     OPEN(UNIT=1,FILE="DATA/MODIS/IND_1deg_1km_EASE_lat.dat",STATUS='OLD')
      
      DO i=1, col_modis
         READ(1,*) IND_1km_ease_lat(i)
      END DO
      CLOSE(1)

      OPEN(UNIT=1,FILE="DATA/MODIS/MISSING_IND_1km_EASE_lat.dat",STATUS='OLD')
       ALLOCATE(MISSING_IND_1km_ease_lat(1408))
      DO i=1, 1408
         READ(1,*) MISSING_IND_1km_ease_lat(i)
      END DO
      CLOSE(1)

     END SUBROUTINE IND_1deg_1km     
    
END MODULE DATADOMAIN
