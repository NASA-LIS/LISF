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
! !ROUTINE: compute_stnwts
! \label{compute_stnwts}
! 
! !REVISION HISTORY: 
!   07-15-05 Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine compute_stnwts(nstns, gridDesc,stnlat, stnlon,&
     npts, stnwt)
  use LIS_logMod, only : LIS_endrun
  implicit none
! !ARGUMENTS:   
  integer     :: nstns
  real        :: gridDesc(50)
  real        :: stnlat(nstns)
  real        :: stnlon(nstns)
  integer     :: npts
  real        :: stnwt(npts,nstns)
!
! !DESCRIPTION:
! 
!  This routine compute the interpolation weights to be applied to 
!  a network of stations, to generate a gridded field from observations
!  from stations. It simply uses an inverse distance weighting 
!  algorithm to compute the relative weights of each station. 
! 
!  \begin{description}
!    \item[nstns]
!     number of stations used in interpolation
!    \item[gridDesc]
!     grid description parameters 
!    \item[stnlat]    
!     station latitudes in degrees 
!    \item[stnlon]    
!     station longitudes in degrees
!    \item[npts]
!     integer maximum number of coordinates
!    \item[stnwt]
!     interpolation weights of stations with respect to the 
!     each point
!    \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[compute\_earth\_coord](\ref{compute_earth_coord}) \newline
!     computes the earth coordinates 
!  \end{description}
!
!EOP

  real :: xpts(npts), ypts(npts)
  real :: rlat(npts), rlon(npts)
  real, parameter     :: fill = -9999.0
  integer :: nv
  integer :: i, j
  real :: dist

  stnwt = 0.0 
!  if(gridDesc(1) .eq. 0) then 
     call compute_earth_coord(gridDesc,npts,fill,xpts,ypts,rlon,rlat,nv,.true.)
     ! compute_earth_coord returns rlon in [0,360].
     ! stnlon is in [-180,180]
     ! This causes a problem when computing the distance.
     ! Translate rlon to [-180,180].
     where ( rlon > 180 ) rlon = rlon - 360
     do i=1,npts
        do j=1,nstns
!           print*, i, rlat(i), stnlat(j), rlon(i), stnlon(j)
           dist = sqrt((rlat(i)-stnlat(j))**2+&
                (rlon(i)-stnlon(j))**2)
           if(dist.eq. 0) then 
              stnwt(i,j) = 1.0
           else
              stnwt(i,j) = 1/(dist)
           endif
        enddo
     enddo
!  else
!     print*, "Sorry, Support for this projection is not implemented yet!"
!     call LIS_endrun
!  endif

end subroutine compute_stnwts
