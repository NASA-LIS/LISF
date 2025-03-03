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
! !ROUTINE:
!  \label{compute_earth_coord_latlon}
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
subroutine compute_earth_coord_latlon(gridDesc,npts,fill,xpts,ypts,& 
     rlon,rlat,nret)

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
!  the specified domain for an equidistant cylindrical projection.
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
  integer :: nret
  real :: rlat1,rlon1,rlat2,rlon2,dlon,dlat
  real :: xmin,xmax,ymin,ymax
  integer :: im,jm,iret,n
  if(gridDesc(1).eq.0) then
     im=gridDesc(2)
     jm=gridDesc(3)
     rlat1=gridDesc(4)
     rlon1=gridDesc(5)
     rlat2=gridDesc(7)
     rlon2=gridDesc(8)
     if(rlat1.gt.rlat2) then 
        dlat=-gridDesc(10)
     else
        dlat=gridDesc(10)
     endif
     if(rlon1.gt.rlon2) then 
        dlon=-gridDesc(9)
     else
        dlon = gridDesc(9)
     endif
     xmin=0
     xmax=im+1
     if(im.eq.nint(360/abs(dlon))) xmax=im+2
     ymin=0
     ymax=jm+1
     nret=0

!  translate grid coordinates to earth coordinates
     do n=1,npts
        if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. & 
             ypts(n).ge.ymin.and.ypts(n).le.ymax) then
           rlon(n)=rlon1+dlon*(xpts(n)-1)
           if(rlon(n).lt.0) then 
              rlon(n) = 360+rlon(n)
           endif
           rlat(n)=rlat1+dlat*(ypts(n)-1)
           nret=nret+1
        else
           rlon(n)=fill
           rlat(n)=fill
        endif
     enddo


!  projection unrecognized
  else
     iret=-1
     do n=1,npts
        rlon(n)=fill
        rlat(n)=fill
     enddo
  endif
end subroutine compute_earth_coord_latlon
