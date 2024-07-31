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
! !ROUTINE: compute_earth_coord_merc
! \label{compute_earth_coord_merc}
!
! !REVISION HISTORY: 
!   04-10-96 Mark Iredell;  Initial Specification
!   07-15-05 Sujay Kumar; Modified verision with floating point arithmetic. 
!
! !INTERFACE:
subroutine compute_earth_coord_merc(gridDesc,npts,fill,xpts,ypts,& 
           rlon,rlat,nret)
! !USES:   
  use map_utils
  implicit none

! !ARGUMENTS: 
  real            :: gridDesc(20)
  integer         :: npts
  real            :: fill
  real            :: xpts(npts),ypts(npts)
  real            :: rlat(npts)
  real            :: rlon(npts)
  integer         :: nret

! !DESCRIPTION:
!  This subroutine computes the earth coordinates of 
!  the specified domain for a mercator projection
!  This routine is based on the
!  decoding routines in the NCEP interoplation package and 
!  has been modified the adopted module from the Weather
!  Research and Forecasting (WRF) model. 
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
!    \item[nret]
!     return code (0-success)
!    \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[map\_set](\ref{map_set}) \newline
!     Sets the projection to mercator
!   \item[ij\_to\_latlon](\ref{ij_to_latlon}) \newline
!     Computes the lat lon values for each i,j
!  \end{description}
!EOP

  type(proj_info) :: proj
  integer :: i

  if(griddesc(1).eq.1) then
     call map_set(PROJ_MERC,gridDesc(4),gridDesc(5),&
          gridDesc(8)*1000,gridDesc(11),gridDesc(10),0.0,&
          nint(gridDesc(2)),nint(gridDesc(3)),proj)
     do i=1,npts
        call ij_to_latlon(proj,xpts(i),ypts(i),rlat(i),rlon(i))
     enddo
  endif
end subroutine compute_earth_coord_merc
