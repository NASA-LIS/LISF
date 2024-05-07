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
! !ROUTINE: compute_earth_coord_latlon_ur
!  \label{compute_earth_coord_latlon_ur}
!
! !REVISION HISTORY: 
!   20 Sept 08 Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine compute_earth_coord_latlon_ur(griddesco,gridDesci,npts,fill,xpts,ypts,& 
     rlon,rlat,nret)

  implicit none

! !ARGUMENTS: 
  real            :: griddesco(50)
  real            :: gridDesci(50)
  integer         :: npts
  real            :: fill
  real            :: xpts(npts),ypts(npts)
  real            :: rlat(npts)
  real            :: rlon(npts)

! !DESCRIPTION:
!  This subroutine computes the upper right earth coordinates of 
!  the specified domain for an equidistant cylindrical projection.
!  This routine is based on the grid
!  decoding routines in the NCEP interoplation package. 
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
  real  :: tmp_lon
  real  :: lat_ll, lon_ll, lat_ur, lon_ur

  lat_ll = gridDesci(4)
  lon_ll = gridDesci(5)
  lat_ur = gridDesci(7)
  lon_ur = gridDesci(8)

!  if(lon_ll.lt.0) then 
!     lon_ll = lon_ll + 360
!  endif
!  if(lon_ur.lt.0) then 
!     lon_ur = lon_ur + 360
!  endif
  if(gridDesco(1).eq.0) then
     im=gridDesco(2)
     jm=gridDesco(3)
     rlat1=gridDesco(4)
     rlon1=gridDesco(5)
     rlat2=gridDesco(7)
     rlon2=gridDesco(8)
     if(rlat1.gt.rlat2) then 
        dlat=-gridDesco(9)
     else
        dlat=gridDesco(9)
     endif
     if(rlon1.gt.rlon2) then 
        dlon=-gridDesco(10)
     else
        dlon = gridDesco(10)
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
           if(rlon(n).gt.180) then 
              tmp_lon = rlon(n)-360
           else
              tmp_lon = rlon(n)
           endif
           rlat(n) = min(rlat(n) + dlat/2,lat_ur)
!           rlon(n) = min(rlon(n) + dlon/2,lon_ur)
           rlon(n) = min(tmp_lon+ dlon/2,lon_ur)
           if(tmp_lon.lt.0) then 
              rlon(n) = 360+tmp_lon
           else
              rlon(n) = tmp_lon
           endif
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
end subroutine compute_earth_coord_latlon_ur
