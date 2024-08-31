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
! !ROUTINE: ezlh_inverse
! \label{ezlh_inverse}
!
! !INTERFACE: 
integer function ezlh_inverse (grid, r, s, lat, lon)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
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
! 
!EOP
!BOP
! 
  implicit none
! !ARGUMENTS: 
  character*(*) grid
  real r, s, lat, lon

! !DESCRIPTION: 
!       convert azimuthal equal area or equal area cylindrical 
!       grid coordinates to geographic coordinates (spherical earth)
!
!       \begin{verbatim}
!       status = ezlh_inverse (grid, r, s, lat, lon)
!
!       input : grid - projection name '[NSM][lh]'
!               where l = "low"  = 25km resolution
!                     h = "high" = 12.5km resolution
!               r, s - grid column and row coordinates
!
!       output: lat, lon - geo. coords. (decimal degrees)
!
!       result: status = 0 indicates normal successful completion
!                       -1 indicates error status (point not on grid)
!
!       \end{verbatim}
!EOP
  integer cols, rows, scale1
  real Rg, phi, lam, rho
  real gamma, beta, epsilon1, x, y
  real sinphi1, cosphi1
  real r0, s0, c
! radius of the earth (km), authalic sphere based on International datum 
  real, parameter :: RE_km = 6371.228
! nominal cell size in kilometers
  real, parameter :: CELL_km = 25.067525

! scale factor for standard paralles at +/-30.00 degrees
  real, parameter :: COS_PHI1 = .866025403

  real, parameter :: PI = 3.141592653589793

  ezlh_inverse = -1
  
  if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then
     cols = 721
     rows = 721
  else if (grid(1:1).eq.'M') then
     cols = 1383
     rows = 586
  else
     print *, 'ezlh_inverse: unknown projection: ', grid
     return
  endif
  
  if (grid(2:2).eq.'l') then
     scale1 = 1
  else if (grid(2:2).eq.'h') then
     scale1 = 2
  else
     print *, 'ezlh_inverse: unknown projection: ', grid
     return
  endif
  
  Rg = scale1 * RE_km/CELL_km
  
  r0 = (cols-1)/2. * scale1
  s0 = (rows-1)/2. * scale1
  
  x = r - r0
  y = -(s - s0)
  
  if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then 
     rho = sqrt(x*x + y*y)
     if (rho.eq.0.0) then
        if (grid(1:1).eq.'N') lat = 90.0 
        if (grid(1:1).eq.'S') lat = -90.0 
        lon = 0.0
     else
        if (grid(1:1).eq.'N') then
           sinphi1 = sin(PI/2.)
           cosphi1 = cos(PI/2.)
           if (y.eq.0.) then
              if (r.le.r0) lam = -PI/2.
              if (r.gt.r0) lam = PI/2.
           else
              lam = atan2(x,-y)
           endif
        else if (grid(1:1).eq.'S') then
           sinphi1 = sin(-PI/2.)
           cosphi1 = cos(-PI/2.)
           if (y.eq.0.) then
              if (r.le.r0) lam = -PI/2.
              if (r.gt.r0) lam = PI/2.
           else
              lam = atan2(x,y)
           endif
        endif
        gamma = rho/(2 * Rg)
        if (abs(gamma).gt.1.) return
        c = 2 * asin(gamma)
        beta = cos(c) * sinphi1 + y * sin(c) * (cosphi1/rho)
        if (abs(beta).gt.1.) return
        phi = asin(beta)
        lat = phi*180.0/PI
        lon = lam*180.0/PI
     endif
     
  else if (grid(1:1).eq.'M') then
     !
     !    allow .5 cell tolerance in arcsin function
     !    so that grid coordinates which are less than .5 cells
     !    above 90.00N or below 90.00S are given a lat of 90.00
     !
     epsilon1 = 1 + 0.5/Rg
     beta = y*COS_PHI1/Rg
     if (abs(beta).gt.epsilon1) return
     if (beta.le.-1.) then
        phi = -PI/2.
     else if (beta.ge.1.) then
        phi = PI/2.
     else
        phi = asin(beta)
     endif
     lam = x/COS_PHI1/Rg
     lat = phi*180.0/PI
     lon = lam*180.0/PI
  endif
  
  ezlh_inverse = 0
  return
end function ezlh_inverse



