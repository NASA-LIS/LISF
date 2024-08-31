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
! !ROUTINE: latlonTopolar
! \label{latlonTopolar}
!
! !REVISION HISTORY: 
!   86-07-17  MCDONELL,J.
!   88-06-07  R.E.JONES   CLEAN UP CODE, TAKE OUT GOTO, USE THEN, ELSE
!   89-11-02  R.E.JONES   CHANGE TO CRAY CFT77 FORTRAN
!   05-27-04  Sujay Kumar Incorporated in LIS
! 
! !INTERFACE:
subroutine latlontopolar(alat,along,xmeshl,orient,xi,xj)

  implicit none
! !ARGUMENTS: 
  real     :: alat
  real     :: along  
  real     :: xmeshl
  real     :: orient
  real     :: xi
  real     :: xj

!
! !DESCRIPTION:
!   Converts the coordinates of a location on earth from the
!   natural coordinate system of latitude/longitude to the grid (i,j)
!   coordinate system overlaid on a polar stereographic map pro-
!   jection true at 60 degrees n or s latitude. 
!
!  The arguments are: 
!  \begin{description}
!  \item[alat]
!    latitude in degrees
!  \item[along]
!    longitude in degrees
!  \item[xmeshl]
!    mesh length of grid in km at 60deg lat (<0 if sh)
!  \item[orient]
!    orientation west longitude of the grid 
!   \item[xi]
!    I of the point relative to North or South pole
!   \item[xj]
!    J of the point relative to North or South pole
!  \end{description}
!
!EOP

  real :: re,xlat,wlong,r
  real,parameter :: radpd = .01745329
  real,parameter :: earthr = 6371.2

  
  RE    = (EARTHR * 1.86603) / XMESHL
  XLAT  = ALAT * RADPD
  
  IF (XMESHL.GE.0.) THEN
     WLONG = (ALONG + 180.0 - ORIENT) * RADPD
     R     = (RE * COS(XLAT)) / (1.0 + SIN(XLAT))
     XI    =   R * SIN(WLONG)
     XJ    =   R * COS(WLONG)
  ELSE
     RE    = -RE
     XLAT  = -XLAT
     WLONG = (ALONG - ORIENT) * RADPD
     R     = (RE * COS(XLAT)) / (1.0 + SIN(XLAT))
     XI    =   R * SIN(WLONG)
     XJ    =  -R * COS(WLONG)
  ENDIF
  
  RETURN
END subroutine latlonTopolar
