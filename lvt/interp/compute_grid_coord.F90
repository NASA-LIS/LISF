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
! !ROUTINE: compute_grid_coord
!  \label{compute_grid_coord}
!
! !INTERFACE:
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
!   04-10-96 Mark Iredell;  Initial Specification
!   05-27-04 Sujay Kumar; Modified verision with floating point arithmetic. 
! 
!EOP
!BOP
!
! !INTERFACE:
subroutine compute_grid_coord(gridDesc,npts,fill,xpts,ypts,rlon,rlat,nret)
  implicit none
! !ARGUMENTS: 
  real        :: gridDesc(50)
  integer     :: npts
  real        :: fill
  real        :: xpts(npts),ypts(npts)
  real        :: rlat(npts)
  real        :: rlon(npts)
  integer     :: nret
! !DESCRIPTION: 
!  This subroutine computes the grid coordinates (cartesian) of 
!  the specified domain. This routine is based on the grid
!  decoding routines in the ipolates interoplation package. 
!  
!  The input options include :
!  The current code recognizes the following projections:
!             (gridDesc(1)=000) equidistant cylindrical
!             (gridDesc(1)=001) mercator cylindrical
!             (gridDesc(1)=003) lambert conformal conical
!             (gridDesc(1)=004) gaussian cylindrical
!             (gridDesc(1)=005) polar stereographic azimuthal
!
!  \begin{description}
!    \item[gridDesc]
!     grid description parameters 
!    \item[npts]
!     integer maximum number of coordinates
!    \item[fill]
!     fill value to set invalid output data
!    \item[xpts]
!     output grid x point coordinates
!    \item[ypts]
!     output grid y point coordinates
!    \item[rlat]    
!     input latitudes in degrees
!    \item[rlon]    
!     input longitudes in degrees
!    \item[nret]
!     return code (0-success)
!    \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[compute\_grid\_coord\_latlon](\ref{compute_grid_coord_latlon} \newline
!     computes the grid coordinates of a latlon grid
!   \item[compute\_grid\_coord\_merc](\ref{compute_grid_coord_merc} \newline
!     computes the grid coordinates of a mercator grid
!   \item[compute\_grid\_coord\_lambert](\ref{compute_grid_coord_lambert} \newline
!     computes the grid coordinates of a lambert conformal grid
!   \item[compute\_grid\_coord\_gauss](\ref{compute_grid_coord_gauss} \newline
!     computes the grid coordinates of a gaussian cylindrical grid
!   \item[compute\_grid\_coord\_polar](\ref{compute_grid_coord_polar} \newline
!     computes the grid coordinates of a polar stereographic grid
!  \end{description}
!EOP

!  equidistant cylindrical
  if(gridDesc(1).eq.0) then
     call compute_grid_coord_latlon(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!  mercator
  elseif(gridDesc(1).eq.1) then      
     call compute_grid_coord_merc(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!  Lambert Conformal
  elseif(gridDesc(1).eq.3) then      
     call compute_grid_coord_lambert(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!     gaussian cylindrical
  elseif(gridDesc(1).eq.4) then
     call compute_grid_coord_gauss(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!  Polar Stereographic 
  elseif(gridDesc(1).eq.5) then
     call compute_grid_coord_polar(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!  HRAP
  elseif(gridDesc(1).eq.8) then
     call compute_grid_coord_hrap(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!  EASE cylindrical
  elseif(gridDesc(1).eq.9) then
     call compute_grid_coord_ease(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
  else
     print*, 'Unrecognized Projection .... '
     print*, 'Program stopping ..'
     stop
  endif
end subroutine compute_grid_coord
