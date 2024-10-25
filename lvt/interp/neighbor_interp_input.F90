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
! !ROUTINE: neighbor_interp_input
! \label{neighbor_interp_input}
!
! !INTERFACE:    
subroutine neighbor_interp_input (gridDesci,gridDesco,npts, rlat2,rlon2,n112)
  implicit none
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subprogram performs issues calls to compute the 
!  neighbor information for neighbor search 
!  interpolation,from any grid to any grid for scalar fields. 
!  The grids are defined by their grid description arrays. 
!  
!  The grid description arrays are based on the decoding 
!  schemes used by NCEP. However, in order to remove the integer
!  arithmetic employed in the original ipolates, the routines
!  are rewritten using real number manipulations. The general 
!  structure remains the same. 
!    
!  The current code recognizes the following projections: \newline
!             (gridDesc(1)=0) equidistant cylindrical \newline
!             (gridDesc(1)=1) mercator cylindrical \newline
!             (gridDesc(1)=3) lambert conformal conical \newline
!             (gridDesc(1)=4) gaussian cylindrical (spectral native) \newline
!             (gridDesc(1)=5) polar stereographic azimuthal \newline
!  where gridDesc could be defined for either the input grid or the 
!  output grid. 
!
!  The arguments are: 
!  \begin{description}
!    \item[gridDesci]
!     input grid description parameters 
!    \item[gridDesco]
!     output grid description parameters 
!    \item[npts] 
!     number of points to in the output field
!    \item[rlat2]    
!     output latitudes in degrees
!    \item[rlon2]    
!     output longitudes in degrees
!    \item[n112]
!     index of neighbor points 
!    \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[compute\_earth\_coord](\ref{compute_earth_coord}) \newline
!     Computes the earth coordinates for the output grid
!   \item[compute\_grid\_coord](\ref{compute_grid_coord}) \newline
!     Computes the grid coordinates of the input grid, based
!     on the earth coordinates of the output grid. 
!   \item[get\_field\_pos](\ref{get_field_pos}) \newline
!     computes the field position for a given point
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! !ARGUMENTS:
  real, intent(in) :: gridDesci(50)
  real             :: gridDesco(50)
  integer          :: npts
  real             :: rlat2(npts)
  real             :: rlon2(npts)
  integer          :: n112(npts)
!
!EOP
  integer             :: n
  integer             :: mo, nv 
  real, parameter     :: fill = -9999.0
  real                :: xpts(npts), ypts(npts)

  integer             :: i1, j1
  real                :: xi, yi
  integer             :: get_fieldpos

  mo = npts
  !------------------------------------------------------------------------
  !  Calls the routines to decode the grid description and 
  !  calculates the weights and neighbor information to perform
  !  spatial interpolation. This routine eliminates the need to 
  !  compute these weights repeatedly during interpolation. 
  !------------------------------------------------------------------------
  if(gridDesco(1).ge.0) then
     call compute_earth_coord(gridDesco,mo,fill,xpts,ypts,rlon2,rlat2,nv)
  endif
  call compute_grid_coord(gridDesci,mo,fill,xpts,ypts,rlon2,rlat2,nv)
  
  do n=1,mo
     xi=xpts(n)
     yi=ypts(n)
     if(xi.ne.fill.and.yi.ne.fill) then
        i1=nint(xi)
        j1=nint(yi)
!        i1=xi
!        j1=yi
        n112(n)=get_fieldpos(i1,j1,gridDesci)
!        print*, 'here ',n,xi,yi,i1,j1,n112(n)
     else
        n112(n)=0
     endif
  enddo
  
end subroutine neighbor_interp_input

