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
! !ROUTINE: compute_earth_coord_ease
!  \label{compute_earth_coord_ease}
!
! !REVISION HISTORY: 
!   04-10-96 Mark Iredell;  Initial Specification
!   05-27-04 Sujay Kumar; Modified verision with floating point arithmetic. 
!   06-11-07 Bailing Li; adapted from compute_earth_coord_latlon
!
! !INTERFACE:
subroutine compute_earth_coord_ease(gridDesc,npts,fill,xpts,ypts,& 
     rlon,rlat,nret)

  use map_utils
  use easeV2_utils

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
!  the specified domain for an EASE cylindrical projection.
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
  real :: cornerlat1, cornerlon1, c0,r0
  integer :: im,jm,iret,n
  type(proj_info) :: proj

  if(gridDesc(1).eq.9) then

     im=gridDesc(2)
     jm=gridDesc(3)
     cornerlat1 = gridDesc(4)
     cornerlon1 = gridDesc(5)

     if(gridDesc(9).eq.1) then !'Ml grid'
        print*, 'Ml grid support needs to be implemented'
        stop
     elseif(gridDesc(9).eq.2) then 
        print*, 'Nl grid support needs to be implemented'
        stop
     elseif(gridDesc(9).eq.3) then 
        print*, 'Sl grid support needs to be implemented'
        stop
     elseif(gridDesc(9).eq.4) then !'M36 grid'
        if(gridDesc(11).eq.1) then 
           do n=1,npts
              if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. & 
                   ypts(n).ge.ymin.and.ypts(n).le.ymax) then
                 !get the grid coordinate of the upper left corner of the domain
                 !which is where EASE grid coordinates start 
                 call easeV2_convert('M36',cornerlat1, cornerlon1,c0,r0)

                 call easeV2_inverse('M36',(xpts(n)+c0)*1.0,&
                      (ypts(n)+r0)*1.0,rlat(n),rlon(n))
              endif
           enddo
        else
           call map_set(PROJ_EASEV2,gridDesc(4),gridDesc(5),&
                gridDesc(9),gridDesc(9),gridDesc(9),0.0,&
                nint(gridDesc(2)),nint(gridDesc(3)),proj)
           do n=1,npts
              if(xpts(n).gt.0.and.ypts(n).gt.0) then
                 call ij_to_latlon(proj, xpts(n),ypts(n),rlat(n), rlon(n))
              else
                 rlat(n) = -9999.0
                 rlon(n) = -9999.0
              endif
           enddo
        endif
     elseif(gridDesc(9).eq.5) then !'M09 grid'
        if(gridDesc(11).eq.1) then 
           do n=1,npts
              if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. & 
                   ypts(n).ge.ymin.and.ypts(n).le.ymax) then
                 !get the grid coordinate of the upper left corner of the domain
                 !which is where EASE grid coordinates start 
                 call easeV2_convert('M09',cornerlat1, cornerlon1,c0,r0)
                 
                 call easeV2_inverse('M09',(xpts(n)+c0)*1.0,&
                      (ypts(n)+r0)*1.0,rlat(n),rlon(n))
              endif
           enddo
        else
           call map_set(PROJ_EASEV2,gridDesc(4),gridDesc(5),&
                gridDesc(9),gridDesc(9),gridDesc(9),0.0,&
                nint(gridDesc(2)),nint(gridDesc(3)),proj)
           do n=1,npts
              if(xpts(n).gt.0.and.ypts(n).gt.0) then
                 call ij_to_latlon(proj, xpts(n),ypts(n),rlat(n), rlon(n))
              else
                 rlat(n) = -9999.0
                 rlon(n) = -9999.0
              endif
           enddo
        endif
     endif
  else
     iret=-1
     do n=1,npts
        xpts(n)=fill
        ypts(n)=fill
     enddo     
  endif

!old stuff
#if 0 

  if(gridDesc(1).eq.7) then
     im=gridDesc(2)
     jm=gridDesc(3)
     rlat1=gridDesc(4)
     rlon1=gridDesc(5)
     rlat2=gridDesc(7)
     rlon2=gridDesc(8)
     if(rlat1.gt.rlat2) then 
        dlat=-gridDesc(9)
     else
        dlat=gridDesc(9)
     endif
     if(rlon1.gt.rlon2) then 
        dlon=-gridDesc(10)
     else
        dlon = gridDesc(10)
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
           !get the grid coordinate of the upper left corner of the domain
           !which is where EASE grid coordinates start 
           ul_lat=max(rlat1,rlat2)
           ul_lon=min(rlon1,rlon2)
           call ezlh_convert('Ml',ul_lat,ul_lon,c0,r0)
           call ezlh_inverse('Ml',(xpts(n)+c0)*1.0,(ypts(n)+r0)*1.0,rlat(n),rlon(n))
           if (rlon(n) <0)then
             rlon(n)=360+rlon(n)
           end if
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
#endif
end subroutine compute_earth_coord_ease
