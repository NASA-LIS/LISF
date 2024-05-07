!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: polarToLatLon
! \label{polarToLatLon}
!
! !REVISION HISTORY
!   86-07-17  R.E.JONES
!   89-11-01  R.E.JONES   CHANGE TO CRAY CFT77 FORTRAN
!   05-27-04  Sujay Kumar Incorporated in LIS
!
! !INTERFACE:
subroutine polarToLatLon(xi,xj,xmeshl,orient,alat,along)

  implicit none
! !ARGUMENTS: 
  real    :: xi
  real    :: xj
  real    :: xmeshl
  real    :: orient
  real    :: alat
  real    :: along
! 
! !DESCRIPTION: 
!   Converts the coordinates of a location from the grid(i,j)
!   coordinate system overlaid on the polar stereographic map projec-
!   tion true at 60 degrees N or S latitude to the natural coordinate
!   system of latitude/longitude on the earth.
!
!  The arguments are: 
!  \begin{description}
!   \item[xi]
!    I of the point relative to North or South pole
!   \item[xj]
!    J of the point relative to North or South pole
!  \item[xmeshl]
!    mesh length of grid in km at 60deg lat (<0 if sh)
!  \item[orient]
!    orientation west longitude of the grid 
!  \item[alat]
!    latitude in degrees
!  \item[along]
!    longitude in degrees
!  \end{description}
!EOP

  real    :: gi2,r2
  real :: angle

  real,parameter :: degprd = 57.2957795
  real,parameter :: earthr = 6371.2

  gi2   = ((1.86603 * earthr) / (xmeshl))**2
  r2    = xi * xi + xj * xj

  if (r2.eq.0.0) then
     along = 0.0
     alat  = 90.0
     if (xmeshl.lt.0.0) alat = -alat
     return
  else
     alat  = asin((gi2 - r2) / (gi2 + r2)) * degprd
     angle = degprd * atan2(xj,xi)
     if (angle.lt.0.0) angle = angle + 360.0
  endif

  if (xmeshl.ge.0.0) then
     along = 270.0 + orient - angle

  else

     along = angle + orient - 270.0
     alat  = -(alat)
  endif

  if (along.lt.0.0)   along = along + 360.0
  if (along.ge.360.0) along = along - 360.0

  return

end subroutine polarToLatLon
