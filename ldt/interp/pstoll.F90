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
! !ROUTINE: pstoll
!  \label{pstoll}
!
! !REVISION HISTORY: 
!     05 aug 1998 initial version....................ssgt mccormick/dnxm
!     10 aug 1999 ported to ibm sp2.  added intent attributes to
!                 arguments................................mr gayno/dnxm
!     25 jul 2005 Adopted in LIS
!   
! !INTERFACE:
subroutine pstoll( hemi, pose, ri, rj, mesh, rlat, rlon )

  implicit none

! !ARGUMENTS: 
  integer,       intent(in)    :: hemi
  integer,       intent(in)    :: pose
  real,          intent(in)    :: ri
  real,          intent(in)    :: rj
  integer,       intent(in)    :: mesh
  real,          intent(out)   :: rlat
  real,          intent(out)   :: rlon

! !DESCRIPTION: 
!     converts polar stereographic grid i/j points to lat/lon
!
!     method: \newline
!     
!     initialize grid specific constants. \newline
!     calculate distance between input point and center point of grid. \newline
!     using this distance, calculate latitude and longitude
!      using trigonometry. \newline
!     adjust sign of longitude according to user preference. \newline
!
!  The arguments are: 
!  \begin{description}
!    \item[hemi]
!     index of the hemisphere (1-nh, 2-sh)
!    \item[pose]
!      longitude increment orientation flag (1-positive east, 
!      0-positive west)
!    \item[ri]
!     i-coordinate 
!    \item[rj]
!     j-coordinate
!    \item[mesh]
!     mesh factor
!    \item[rlat]    
!     output latitude in degrees
!    \item[rlon]    
!     output longitude in degrees
!    \end{description}
!EOP
  real,          parameter     :: dierth = 6371.2213
  real,          parameter     :: rad2dg = 57.29577951308
  real,          parameter     :: sclprm = 1.8660254037844
  real,          parameter     :: d2s2   = (dierth * dierth) * &
       (sclprm * sclprm) 
  
  integer                      :: numpts

  
  real                         :: cval
  real                         :: grddis
  real                         :: h(2)
  real                         :: latdeg
  real                         :: latrad
  real                         :: midpt
  real                         :: x
  real                         :: x2y2
  real                         :: y
  
  data h  / 1.0, -1.0 /
!-------------------------------------------------------------------
! the grid distance is determined by dividing the grid distance of a
! whole-mesh grid by the mesh factor.
!-------------------------------------------------------------------

  grddis = 381.0 / real( mesh )
  cval   = d2s2 / (grddis ** 2)

!-------------------------------------------------------------------
! calculate distance from center point
!-------------------------------------------------------------------

  numpts = (64 * mesh) + 1
  midpt  = real( (numpts / 2) + 1 )
  x      = ri - midpt
  y      = -h(hemi) * (rj - midpt)

!-------------------------------------------------------------------
! calculate (x**2 + y**2)
!-------------------------------------------------------------------

  x2y2 = ((x * x) + (y * y))

!-------------------------------------------------------------------
! calculate the latitude angle
!-------------------------------------------------------------------
  
  latrad = (cval - x2y2) / (cval + x2y2)

!-------------------------------------------------------------------
! convert the latitude angle to degrees
!-------------------------------------------------------------------

  latdeg = asin( latrad ) * rad2dg

!-------------------------------------------------------------------
! the latitude of the angle is negative if the angle is in the
! southern hemisphere (hemi = 2)
!-------------------------------------------------------------------

  rlat = latdeg * h(hemi)

! -------------------------------------------------------------------
! calculate the longitude angle of the grid point in radians, then
! convert to degrees
!-------------------------------------------------------------------

  if( y .gt. 0.0 )then

!-------------------------------------------------------------------
! the point is in the western hemisphere (+ 10 degrees)
!-------------------------------------------------------------------

     rlon = 10.0 + ( acos( x / sqrt(x2y2) ) * rad2dg )

  elseif( y .lt. 0.0 )then

!-------------------------------------------------------------------
! the point is in the eastern hemisphere (+ 10 degrees)
!-------------------------------------------------------------------

     rlon = 10.0 - ( acos( x / sqrt(x2y2) ) * rad2dg )

  else

! -------------------------------------------------------------------
! the point is along the +10 or +190 degree line
! -------------------------------------------------------------------

     if( x .gt. 0 )then

! -------------------------------------------------------------------
! the point is on the +10 degree line
!-------------------------------------------------------------------

        rlon = 10.0

     elseif( x .lt. 0 )then

! -------------------------------------------------------------------
!   the point is on the +190.0 degree line
!-------------------------------------------------------------------

        rlon = 190.0
        
     else
           
!-------------------------------------------------------------------
! the point is on the pole, so set the longitude to 0.0
!-------------------------------------------------------------------

        rlon = 0.0

     endif
  endif
      
!-------------------------------------------------------------------
! make all longitude values positive
!-------------------------------------------------------------------

  if( rlon .lt. 0.0 )then
     rlon = 360.0 + rlon
  endif

!-------------------------------------------------------------------
! adjust for positive west (if applicable)
!-------------------------------------------------------------------

  if( pose .ne. 1 )then
     rlon = 360.0 - rlon
  endif
  
  return
end subroutine pstoll
