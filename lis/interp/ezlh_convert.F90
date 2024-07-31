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
! !ROUTINE: ezlh_convert
! \label{ezlh_convert}
! 
! !REVISION HISTORY: 
! 
!       30-Jan.-1992 H.Maybee
!       20-Mar-1992 Ken Knowles  303-492-0644  knowles@kryos.colorado.edu
!       16-Dec-1993 MJ Brodzik   303-492-8263  brodzik@jokull.colorado.edu
!                   Copied from nsmconv.f, changed resolutions from 
!                   40-20-10 km to 25-12.5 km
!       21-Dec-1993 MJ Brodzik   303-492-8263  brodzik@jokull.colorado.edu
!                   Fixed sign of Southern latitudes in ease_inverse.
!       12-Sep-1994 David Hoogstrate 303-492-4116 hoogstra@jokull.colorado.edu
!                   Changed grid cell size. Changed "c","f" to "l","h"
!       25-Oct-1994 David Hoogstrate 303-492-4116 hoogstra@jokull.colorado.edu
!                   Changed row size from 587 to 586 for Mercator projection
!                   Changed function names to "ezlh-.."
!       07-11-2007  Bailing Li; downloaded this code from NSIDC 
! !INTERFACE: 
integer function ezlh_convert (grid, lat, lon, r, s)

  implicit none
! !ARGUMENTS: 
  character*(*) grid
  real lat, lon, r, s

! !DESCRIPTION: 
!       convert geographic coordinates (spherical earth) to 
!       azimuthal equal area or equal area cylindrical grid coordinates
!
!       \begin{verbatim}
!       status = ezlh_convert (grid, lat, lon, r, s)
!
!       input : grid - projection name '[NSM][lh]'
!               where l = "low"  = 25km resolution
!                     h = "high" = 12.5km resolution
!               lat, lon - geo. coords. (decimal degrees)
!
!       output: r, s - column, row coordinates
!
!       result: status = 0 indicates normal successful completion
!                       -1 indicates error status (point not on grid)
!       \end{verbatim}
!EOP

  integer cols, rows, scale
  real Rg, phi, lam, rho
  real r0, s0
!       radius of the earth (km), authalic sphere based on International datum 
  real, parameter :: RE_km = 6371.228
!       nominal cell size in kilometers
  real, parameter :: CELL_km = 25.067525

!       scale factor for standard paralles at +/-30.00 degrees
  real, parameter :: COS_PHI1 = .866025403

  real, parameter :: PI = 3.141592653589793
        
  ezlh_convert = -1
  
  if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then
     cols = 721
     rows = 721
  else if (grid(1:1).eq.'M') then
     cols = 1383
     rows = 586
  else
     print *, 'ezlh_convert: unknown projection: ', grid
     return
  endif

  if (grid(2:2).eq.'l') then
     scale = 1
  else if (grid(2:2).eq.'h') then
     scale = 2
  else
     print *, 'ezlh_convert: unknown projection: ', grid
     return
  endif

  Rg = scale * RE_km/CELL_km

!
!       r0,s0 are defined such that cells at all scales
!       have coincident center points
!
  r0 = (cols-1)/2. * scale
  s0 = (rows-1)/2. * scale
  
  phi = lat*PI/180.
  lam = lon*PI/180.
  
  if (grid(1:1).eq.'N') then
     rho = 2 * Rg * sin(PI/4. - phi/2.)
     r = r0 + rho * sin(lam)
     s = s0 + rho * cos(lam)
     
  else if (grid(1:1).eq.'S') then
     rho = 2 * Rg * cos(PI/4. - phi/2.)
     r = r0 + rho * sin(lam)
     s = s0 - rho * cos(lam)
     
  else if (grid(1:1).eq.'M') then
     r = r0 + Rg * lam * COS_PHI1
     s = s0 - Rg * sin(phi) / COS_PHI1
     
  endif
  
  ezlh_convert = 0
  return
end function ezlh_convert



