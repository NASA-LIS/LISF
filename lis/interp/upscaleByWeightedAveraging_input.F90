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
! !ROUTINE: upscaleByWeightedAveraging_input
! \label{upscaleByWeightedAveraging_input}
!
! !INTERFACE:    
subroutine upscaleByWeightedAveraging_input (gridDesci,gridDesco,mi,mo,& 
     n11,w11)

  implicit none 
! !ARGUMENTS:
  real, intent(in)    :: gridDesci(50)
  real                :: gridDesco(50)
  integer             :: mi
  integer             :: mo
  integer             :: n11(mi)
  real                :: w11(mi)
!
! !DESCRIPTION: 
!  This subprogram performs issues calls to compute the 
!  neighbor information for upscaling data for scalar fields
!  The grids are defined by their grid description arrays. 
!  
!  The basic logic is as follows: The code first computes the 
!  earth coordinates (lat/lon) for the input grid (fine resolution
!  grid). Then it looks up the grid coordinates (x,y) of each 
!  lat/lon in the output grid (coarse resolution grid). The code
!  stores these mappings in the n11 array, which represents the 
!  location of each input grid point in the output grid. 
!
!  Note that the input grid needs to be setup either as a global 
!  grid or with halos padded around it (if using a subset of the
!  global grid), for the algorithm to work properly in a 
!  multiprocessor environment.This is to avoid a case where the 
!  required input grid points are split across different processors.  
!  
!  The grid description arrays are based on the decoding 
!  schemes used by NCEP. 
!
!  The current code recognizes the following projections: \newline
!             (gridDesc(1)=0) equidistant cylindrical \newline
!  where gridDesc could be defined for either the input grid or the 
!  output grid. 
!
!  The arguments are: 
!  \begin{description}
!    \item[gridDesci]
!     input grid description parameters 
!    \item[gridDesco]
!     output grid description parameters 
!    \item[mi] 
!     total number of points in the input grid
!    \item[mo] 
!     total number of points in the output grid
!    \item[n11] 
!     array that maps the location of each input grid
!     point in the output grid. 
!    \end{description}
!
!
!EOP

  real                :: xcompact, ycompact
  integer             :: n 
  integer             :: i,j
  real                :: xi, yi
  real                :: xpts(mi), ypts(mi)
  real                :: rlat(mi), rlon(mi)

  integer             :: nv
  real                :: dx,dy,dist,maxdist
  integer             :: get_fieldpos
  real, parameter     :: fill = -9999.0

  w11 = -9999.0
  xcompact = 0.25
  ycompact = 0.25

  maxdist = xcompact**2 + ycompact**2

  if(gridDesci(1).ge.0) then 
! Compute the lat/lons on the fine resolution grid. 
     call compute_earth_coord(gridDesci,mi,fill,xpts,ypts,rlon,rlat,nv,.true.)
  endif
! Compute the xy values on the coarse resolution grid
  call compute_grid_coord(gridDesco,mi,fill, xpts, ypts, rlon, rlat, nv)

!loop through the input grid points
  do n=1,mi
     xi = xpts(n)
     yi = ypts(n)
     dx = xi -nint(xi)
     dy = yi -nint(yi)
     
     if(xi.ne.fill.and.yi.ne.fill) then 
        i = nint(xi)
        j = nint(yi)
        n11(n) = get_fieldpos(i,j,gridDesco)

        dist   = dx**2 + dy**2
        
        if(dist.le.maxdist) then 
           w11(n) = (maxdist - dist)/(maxdist+dist)
        else
           w11(n) = 0.0
        endif
     else
        n11(n) = 0
        w11(n) = 0 
     endif
  enddo

end subroutine upscaleByWeightedAveraging_input


