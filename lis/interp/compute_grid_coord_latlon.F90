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
! !ROUTINE: compute_grid_coord_latlon
!  \label{compute_grid_coord_latlon}
!
! !REVISION HISTORY: 
!   04-10-96 Mark Iredell;  Initial Specification
!   05-27-04 Sujay Kumar; Modified verision with floating point arithmetic. 
!
! !INTERFACE:
subroutine compute_grid_coord_latlon(gridDesc,npts,fill,xpts,ypts,& 
     rlon,rlat,nret)

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
!  the specified domain for an equidistant cylindrical rojection.
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
!EOP
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

     do n=1,npts
        if(abs(rlon(n)).le.360.and.abs(rlat(n)).le.90) then
           if(rlon(n).gt.180) then 
              xpts(n)=1+(rlon(n)-360-rlon1)/dlon
           else
              xpts(n) = 1+(rlon(n)-rlon1)/dlon
           endif
           ypts(n)=1+(rlat(n)-rlat1)/dlat

           if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. & 
                ypts(n).ge.ymin.and.ypts(n).le.ymax) then
              nret=nret+1
           else
              xpts(n)=fill
              ypts(n)=fill
           endif
        else
           xpts(n)=fill
           ypts(n)=fill
        endif
     enddo
  else
     iret=-1
  endif
end subroutine compute_grid_coord_latlon
