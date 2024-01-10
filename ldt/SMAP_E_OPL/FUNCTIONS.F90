!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
MODULE FUNCTIONS

CONTAINS

    FUNCTION LAT(lat_up,lat_lo,space) RESULT(lat_vec)
    
      IMPLICIT NONE
      INTEGER       :: i, lat_count
      REAL*8 :: lat_up, lat_lo, space
      REAL*8,DIMENSION(:),ALLOCATABLE:: lat_vec      
      lat_count = (lat_up-lat_lo)/space 
      ALLOCATE(lat_vec(lat_count))
      DO i=1,lat_count
         lat_vec(i) = lat_up-space/2.0-(i-1)*space
      END DO
    END FUNCTION LAT

    FUNCTION LON(lon_lf,lon_rt,space) RESULT(lon_vec)

      IMPLICIT NONE
      INTEGER       :: i, lon_count
      REAL*8 :: lon_lf,lon_rt,space
      REAL*8,DIMENSION(:),ALLOCATABLE:: lon_vec

      lon_count = (lon_rt-lon_lf)/space
      ALLOCATE(lon_vec(lon_count))
      DO i=1,lon_count
         lon_vec(i) = lon_lf+space/2.0+(i-1)*space
      END DO

    END FUNCTION LON

   FUNCTION SEQ(start,endnum, interval) RESULT(seq_vec)
      IMPLICIT NONE
      INTEGER       :: i, num_count
      REAL*4 :: start, endnum, interval
      REAL*4,DIMENSION(:),ALLOCATABLE:: seq_vec
      num_count = (start-endnum)/interval+1
      ALLOCATE(seq_vec(num_count))
      DO i=1,num_count
         seq_vec(i) = start+(i-1)*interval
      END DO
   END FUNCTION SEQ

  FUNCTION DECYEAR(year,doy) RESULT(yeardec)
      IMPLICIT NONE
      INTEGER    :: leap, year, doy
      REAL*8     :: yeardec
      leap=MOD(year,4)
      IF (leap.eq.0) THEN
         !leap year
         yeardec=REAL(year,8)+REAL(doy-1,8)/366.0
      
      ELSE
         yeardec=REAL(year,8)+REAL(doy-1,8)/365.0
      END IF
  END FUNCTION DECYEAR
 
  FUNCTION MONOFYR(year,doy) RESULT(mon)
  IMPLICIT NONE
      INTEGER    :: leap, year, doy, mon, m, days
      INTEGER    :: month(1:12)
      
      leap=MOD(year,4)

      IF (leap.eq.0) THEN
          !leap year
          month(1:12) = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
         
      ELSE 
          month(1:12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      END IF
      days=0 
      DO m=1,12
         days=month(m)+days
         IF (doy.le.days) THEN
            mon=m
            EXIT
         END IF
       !PRINT *,'doy days m mom', doy, days, m, mon
      END DO
  END FUNCTION MONOFYR


  FUNCTION DAYOFYR(year,mon,day) RESULT(doy)
  IMPLICIT NONE
      INTEGER    :: leap, year, day, doy, mon, m, sumdays
      INTEGER    :: month(1:12)

      leap=MOD(year,4)

      IF (leap.eq.0) THEN
          !leap year
          month(1:12) = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

      ELSE
          month(1:12) = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      END IF
      sumdays=0
      IF (mon.EQ.1) THEN
          doy=day
      ELSE
           DO m=2,mon
             sumdays=sumdays+month(m-1)
           END DO
           doy=sumdays+day
      END IF
  END FUNCTION DAYOFYR

END MODULE FUNCTIONS
