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
! !ROUTINE: compute_earth_coord_UTM
!  \label{compute_earth_coord_UTM}
!
! !REVISION HISTORY: 
!   04-10-96 Mark Iredell;  Initial Specification
!   05-27-04 Sujay Kumar; Modified verision with floating point arithmetic. 
!
! !INTERFACE:
subroutine compute_earth_coord_UTM(gridDesc,npts,fill,xpts,ypts,& 
     rlon,rlat,nret)

  use UTM_utils,   only : UTM2geo

  implicit none

! !ARGUMENTS: 
  real            :: gridDesc(50)
  integer         :: npts
  real            :: fill
  real            :: xpts(npts),ypts(npts)
  real            :: rlat(npts)
  real            :: rlon(npts)

! !DESCRIPTION:
!  This subroutine computes the earth coordinates of 
!  the specified domain for a UTM projection
!  
!  \begin{description}
!    \item[gridDesc]
!     grid description parameters 
!    \item[npts]
!     integer maximum number of coordinates
!    \item[fill]
!     fill value to set invalid output data
!    \item[xpts]
!     grid x point coordinates
!    \item[ypts]
!     grid y point coordinates
!    \item[rlat]    
!     output latitudes in degrees
!    \item[rlon]    
!     output longitudes in degrees
!    \end{description}
!
!EOP
  integer              :: nret
  integer              :: n
  real                 :: northing, easting

  do n=1,npts
     northing = gridDesc(4) + (ypts(n)-1)*gridDesc(9)
     easting = gridDesc(5) + (xpts(n)-1)*gridDesc(9)
     call UTM2Geo(nint(gridDesc(10)), northing, easting,rlat(n),rlon(n))
  enddo
end subroutine compute_earth_coord_UTM
