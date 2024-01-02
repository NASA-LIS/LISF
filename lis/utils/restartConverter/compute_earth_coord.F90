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
! !ROUTINE: compute_earth_coord
!  \label{compute_earth_coord}
!
! !REVISION HISTORY: 
!   04-10-96 Mark Iredell;  Initial Specification
!   05-27-04 Sujay Kumar; Modified verision with floating point arithmetic. 
!
! !INTERFACE:
subroutine compute_earth_coord(gridDesc,npts,fill,xpts,ypts,rlon,rlat,nret)

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
!  This subroutine computes the earth coordinates (lat/lon values) 
!  of the specified domain. This routine is based on the grid
!  decoding routines in the ipolates interoplation package. 
!  
!  The input options include :
!  The current code recognizes the following projections: \newline
!             (gridDesc(1)=000) equidistant cylindrical \newline
!             (gridDesc(1)=001) mercator cylindrical \newline
!             (gridDesc(1)=003) lambert conformal conical \newline
!             (gridDesc(1)=004) gaussian cylindrical \newline
!             (gridDesc(1)=005) polar stereographic azimuthal \newline
!
!  \begin{description}
!    \item[gridDesc]
!     grid description parameters 
!    \item[npts]
!     integer maximum number of coordinates
!    \item[fill]
!     fill value to set invalid output data
!    \item[xpts]
!     input grid x point coordinates
!    \item[ypts]
!     input grid y point coordinates
!    \item[rlat]    
!     output latitudes in degrees
!    \item[rlon]    
!     output longitudes in degrees
!    \item[nret]
!     return code (0-success)
!    \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[compute\_earth\_coord\_latlon](\ref{compute_earth_coord_latlon}) \newline
!     computes the earth coordinates of a latlon grid
!   \item[compute\_earth\_coord\_merc](\ref{compute_earth_coord_merc}) \newline
!     computes the earth coordinates of a mercator grid
!   \item[compute\_earth\_coord\_lambert](\ref{compute_earth_coord_lambert}) \newline
!     computes the earth coordinates of a lambert conformal grid
!   \item[compute\_earth\_coord\_gauss](\ref{compute_earth_coord_gauss}) \newline
!     computes the earth coordinates of a gaussian cylindrical grid
!   \item[compute\_earth\_coord\_polar](\ref{compute_earth_coord_polar}) \newline
!     computes the earth coordinates of a polar stereographic grid
!  \end{description}
!EOP
  integer :: im,jm,nm,n
  integer :: i,j

  im=gridDesc(2)
  jm=gridDesc(3)
  nm=im*jm
  if(nm.le.npts) then
     do n=1,nm
        j=(n-1)/im+1
        i=n-im*(j-1)
        xpts(n)=i
        ypts(n)=j
     enddo
     do n=nm+1,npts
        xpts(n)=fill
        ypts(n)=fill
     enddo
  else
     do n=1,npts
        xpts(n)=fill
        ypts(n)=fill
     enddo
  endif

!  equidistant cylindrical
  if(gridDesc(1).eq.0) then
     call compute_earth_coord_latlon(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!  Mercator
  elseif(gridDesc(1).eq.1) then 
     call compute_earth_coord_merc(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!  Lambert Conformal
  elseif(gridDesc(1).eq.3) then 
     call compute_earth_coord_lambert(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!     gaussian cylindrical
  elseif(gridDesc(1).eq.4) then
     call compute_earth_coord_gauss(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!  Polar Stereographic 
  elseif(gridDesc(1).eq.5) then
     call compute_earth_coord_polar(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!  HRAP
  elseif(gridDesc(1).eq.6) then
     call compute_earth_coord_hrap(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
!  EASE cylindrical
  elseif(gridDesc(1).eq.7) then
     call compute_earth_coord_ease(gridDesc,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
  endif
end subroutine compute_earth_coord
