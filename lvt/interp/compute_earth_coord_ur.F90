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
! !ROUTINE: compute_earth_coord_ur
!  \label{compute_earth_coord_ur}
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
!   20 Sept 08 Sujay Kumar; Initial Specification
! 
!EOP
!BOP
!
! !INTERFACE:
subroutine compute_earth_coord_ur(gridDesco,gridDesci,npts,fill,xpts,ypts,rlon,rlat,nret)

  implicit none
! !ARGUMENTS: 
  real        :: gridDesco(50)
  real        :: gridDesci(50)
  integer     :: npts
  real        :: fill
  real        :: xpts(npts),ypts(npts)
  real        :: rlat(npts)
  real        :: rlon(npts)
  integer     :: nret

! !DESCRIPTION: 
!  This subroutine computes the upper right earth coordinates (lat/lon values) 
!  of the specified domain. This routine is based on the grid
!  decoding routines in the ipolates interoplation package. 
!  
!  The input options include :
!  The current code recognizes the following projections: \newline
!             (gridDesc(1)=000) equidistant cylindrical \newline
!
!  \begin{description}
!    \item[gridDesci]
!     input grid description parameters 
!    \item[gridDesco]
!     output grid description parameters 
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
!   \item[compute\_earth\_coord\_latlon\_ur](\ref{compute_earth_coord_latlon_ur}) \newline
!     computes the upper right earth coordinates of a latlon grid
!  \end{description}
!EOP

  integer :: im,jm,nm,n
  integer :: i,j

  im=gridDesco(2)
  jm=gridDesco(3)
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
  if(gridDesco(1).eq.0) then
     call compute_earth_coord_latlon_ur(gridDesco,gridDesci,npts,fill,xpts,ypts,&
          rlon,rlat,nret)
  else
     write(*,*) 'Upper right corner finding algorithm is not implemented'
     write(*,*) 'for this LIS projection'
     write(*,*) 'Program stopping... '
     stop
  endif
end subroutine compute_earth_coord_ur
