!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
MODULE VARIABLES

     IMPLICIT NONE
     INTEGER :: row_ldas, col_ldas, row_avhr, col_avhr, smap_row, smap_col, row_modis, col_modis, smap_geo_row, smap_geo_col
     INTEGER :: row_ease_1km, col_ease_1km, arfs_nrow_lat, arfs_mcol_lon, arfs_nrow_lat_3km, arfs_mcol_lon_3km
     INTEGER,DIMENSION(:),ALLOCATABLE :: IND_1km_ease_lon, IND_1km_ease_lat, MISSING_IND_1km_ease_lat
     REAL*8,DIMENSION(:),ALLOCATABLE :: smap_lat, smap_lon, EASE1KM_lat, EASE1KM_lon
     REAL*4 :: avhr_space, avhr_lat_up, avhr_lat_lo, avhr_lon_lf, avhr_lon_rt, smap_space
     REAL*4 :: ldas_space, ldas_lat_up, ldas_lat_lo, ldas_lon_lf, ldas_lon_rt
     REAL*4 :: modis_space, modis_lat_up, modis_lat_lo, modis_lon_lf, modis_lon_rt
     REAL*4 :: smap_geo_space, smap_geo_lon_lf, smap_geo_lon_rt, smap_geo_lat_up, smap_geo_lat_lo
     REAL*8 :: arfs_geo_lon_lf, arfs_geo_lon_rt, arfs_geo_lat_up, arfs_geo_lat_lo, arfs_lon_space, arfs_lat_space,arfs_lon_3km_space, arfs_lat_3km_space
     REAL*4 :: ROI_lon_lf, ROI_lon_rt, ROI_lat_up, ROI_lat_lo
END MODULE VARIABLES
