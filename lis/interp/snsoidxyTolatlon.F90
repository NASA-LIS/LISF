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
!BOP
!
! !ROUTINE: snsoidxyTolatlon
! \label{snsoidxyTolatlon}
!
! !REVISION HISTORY: 
!   03-14-09  Soni Yatheendradas; Adopted in LIS 5.0 BETA RFS version from GCTP
!
!  CAUTION: This code is specific only to RFS project Sierra Nevada basins zone and not generic 
! !INTERFACE:
subroutine snsoidxyTolatlon(xx,yy,npts,false_easting,false_northing,RRR,lon_center,alat,along)

  implicit none
! !ARGUMENTS: 
  real,intent(in)     :: xx(npts)
  real,intent(in)     :: yy(npts)
  integer,intent(in)  :: npts
  real,intent(in)     :: false_easting
  real,intent(in)     :: false_northing
  real,intent(in)     :: RRR
  real,intent(in)     :: lon_center
  real     :: alat(npts)
  real     :: along(npts) 
 
!
! !DESCRIPTION:
!   Converts the coordinates of a location on earth from the
!   sinusiodal grid (x,y) co-ordinate system of the MODIS tile
!   over the RFS project Sierra Nevada basins, to the natural 
!   coordinate system of latitude/longitude grid 
!   IMPORTANT: This is a VECTORIZED code 
!
!  The arguments are: 
!  \begin{description}
!   \item[xx]
!    xx of the point/s 
!   \item[yy]
!    yy of the point/s 
!   \item[npts]
!    number of point/s 
!   \item[false\_easting]
!    False easting in the same units as the semi-major axis
!   \item[false\_northing]
!    False northing in the same units as the semi-major axis
!   \item[RRR]
!    Radius of the reference sphere
!   \item[lon\_center]
!    Longitude of the central meridian
!  \item[alat]
!    latitude in degrees
!  \item[along]
!    longitude in degrees
!  \end{description}
!
!EOP

  real,PARAMETER     :: PI_HERE=3.141592653589793238
  real,PARAMETER     :: HALF_PI_HERE=PI_HERE*0.5
  real,PARAMETER     :: EPSLN_HERE=1.0e-10
  real               :: xx2(npts)
  real               :: yy2(npts)
  real               :: temptemp(npts)
  integer            :: ii
  
  print*,"Entering snsoidxyTolatlon"
  xx2 = xx - false_easting
  yy2 = yy - false_northing
  alat = yy2 / RRR
  temptemp = abs(alat) - HALF_PI_HERE
  do ii=1,npts
    if (abs(temptemp(ii)) > EPSLN_HERE) then
      along(ii)=lon_center*PI_HERE/180.0+xx2(ii)/(RRR*cos(alat(ii)))
    else   
      along(ii)=lon_center*PI_HERE/180.0
    end if
  end do 
  alat=alat*180.0/PI_HERE
  along=along*180.0/PI_HERE
!  ii=1
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii)
!  ii=235201
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
!  ii=5757601
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
!  ii=124853
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
!  ii=2573378
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
!  ii=599
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
!  ii=23163
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
!  ii=5759553
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
!  ii=5757599
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
!  ii=2400
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
!  ii=2637600
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
!  ii=5760000
!  print*,"tt,xx,yy,alat,along:",ii,xx(ii),yy(ii),alat(ii),along(ii) 
  
  RETURN
END subroutine snsoidxyTolatlon
