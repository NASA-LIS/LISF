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
! !ROUTINE: lltops
!  \label{lltops}
!
! !REVISION HISTORY:
!    05 aug 1998 initial version.....................ssgt mccormick/dnxm
!    10 aug 1999 ported to ibm sp-2.  added intent attributes to
!                arguments.................................mr gayno/dnxm
!
!    29 oct 2005; Sujay Kumar; Incorporated into LIS
!
! !INTERFACE:
subroutine lltops( pose, rlat, rlon, mesh, hemi, ri, rj )

  implicit none

! !ARGUMENTS: 
  integer,      intent(in)     :: pose
  real,         intent(in)     :: rlat
  real,         intent(in)     :: rlon
  integer,      intent(in)     :: mesh
  integer,      intent(out)    :: hemi
  real,         intent(out)    :: ri
  real,         intent(out)    :: rj

! !DESCRIPTION:
!    converts lat/lon to polar stereographic grid i/j points
!
!    method \newline
!
!    calculate map related constants/factors. \newline
!    adjust the longitude according to the user specified
!      sign convention. \newline
!    convert longitude to radians. \newline
!    calculate the distance from the pole to the point in radians. \newline
!    calculate the center point coordinate on the grid. \newline
!    calculate the grid coordinates for the point. \newline
!
!  The arguments are: 
!  \begin{description}
!    \item[pose]
!      longitude increment orientation flag (1-positive east, 
!      0-positive west)
!    \item[rlat]    
!     input latitude in degrees
!    \item[rlon]    
!     input longitude in degrees
!    \item[mesh]
!     mesh factor
!    \item[hemi]
!     index of the hemisphere (1-nh, 2-sh)
!    \item[ri]
!     i-coordinate 
!    \item[rj]
!     j-coordinate
!    \end{description}
!
!EOP

!     ------------------------------------------------------------------
!     set up initial values.  because this routine might be called num-
!     erous times, the following constants have been pre-calculated to
!     save processing time:
!
!         dg2rad = pi / 180.0
!         sclprm = 1.0 + sin(true lat) where trulat is 60 degrees
!         d2s2   = (diam earth ** 2) * (sclprm ** 2)
!                   diam earth is 6371.2213
!     ------------------------------------------------------------------

  real,         parameter      :: dg2rad = 1.745329251994e-2
  real,         parameter      :: dierth = 6371.2213
  real,         parameter      :: sclprm = 1.8660254037844
  real,         parameter      :: d2s2   = (dierth * dierth) * &
       (sclprm * sclprm)

  integer                      :: numpts

  
  real                         :: cval
  real                         :: dist
  real                         :: grddis
  real                         :: latdeg
  real                         :: londeg
  real                         :: lonrad
  real                         :: midpt
  real                         :: polrad


!     ------------------------------------------------------------------
!     the grid distance is determined by dividing the grid distance of a
!     whole-mesh grid by the mesh factor.
!     ------------------------------------------------------------------

  grddis = 381.0 / float(mesh)
  
  cval   = d2s2 / (grddis ** 2)
  
  polrad = sqrt(cval)

!     ------------------------------------------------------------------
!     adjust the longitude increment direction.
!     ------------------------------------------------------------------
  
  if (pose .eq. 1) then
     
     londeg = rlon
     
  else
     
     londeg = 360.0 - rlon
     
  endif

!     ------------------------------------------------------------------
!     calculate longitude in radians.
!     ------------------------------------------------------------------

  if (rlat .ge. 0) then
     
     lonrad = (londeg + 350.0) * dg2rad
     latdeg = rlat
     hemi   = 1

  else
     
     lonrad = (370.0 - londeg) * dg2rad
     latdeg = - rlat
     hemi   = 2
     
  endif

!     ------------------------------------------------------------------
!     calculate the distance from the pole to the point in radians.
!     ------------------------------------------------------------------
  
  dist = polrad * tan((90.0 - latdeg) * dg2rad / 2.0)

!     ------------------------------------------------------------------
!     calculate the center point coordinate on the grid.
!     ------------------------------------------------------------------

  numpts = (64 * mesh) + 1
  
  midpt  = float((numpts / 2) + 1)

!     ------------------------------------------------------------------
!     calculate the grid coordinates for the point.
!     ------------------------------------------------------------------

  ri = midpt + cos(lonrad) * dist
  rj = midpt - sin(lonrad) * dist
  
  return
end subroutine lltops
