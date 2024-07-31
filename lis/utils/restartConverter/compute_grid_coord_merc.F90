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
! !ROUTINE: compute_grid_coord_merc
!  \label{compute_grid_coord_merc}
!
! !REVISION HISTORY: 
!   04-10-96 Mark Iredell;  Initial Specification
!   07-15-05 Sujay Kumar; Modified verision with floating point arithmetic. 
!
! !INTERFACE:
subroutine compute_grid_coord_merc(gridDesc,npts,fill,xpts,ypts,& 
     rlon,rlat,nret)
! !USES:   

  implicit none

! !ARGUMENTS: 
  real            :: gridDesc(50)
  integer         :: npts
  real            :: fill
  real            :: xpts(npts),ypts(npts)
  real            :: rlat(npts)
  real            :: rlon(npts)
  integer         :: nret

! !DESCRIPTION:
!  This subroutine computes the grid coordinates of 
!  the specified domain for a mercator projection.
!  This routine is based on the grid
!  decoding routines in the NCEP interoplation package. 
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
!    \end{description}
!
! !NOTE:
!  This routine is currently unsupported. 
!EOP

  real, parameter:: RERTH=6.3712E3
  real, parameter:: PI=3.14159265358979,DPR=180./PI
  
  write(*,*) &
       'Transformation from mercator projection is not supported'
  stop

end subroutine compute_grid_coord_merc
