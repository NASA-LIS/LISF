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
! !ROUTINE: compute_grid_coord_hrap
!  \label{compute_grid_coord_hrap}
!
! !REVISION HISTORY: 
!   04-10-96 Mark Iredell;  Initial Specification
!   07-15-05 Sujay Kumar; Modified verision with floating point arithmetic. 
!   07-25-06 Sujay Kumar; Updated version to include Map Utils routines
!   02-04-14 Shugong Wang; Fix bugs with LDT_XMRG_Reader 
! 
! !INTERFACE:
subroutine compute_grid_coord_hrap(gridDesc,npts,fill,xpts,ypts,& 
     rlon,rlat,nret)

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
!  This subroutine computes the grid coordinates of 
!  the specified domain for a hrap projection.
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
!  The routines invoked are: 
!  \begin{description}
!   \item[latlonTohrap](\ref{latlonTohrap}) \newline
!     Converts lat, lons to hrap coordinates
!   \item[map\_set](\ref{map_set}) \newline
!     Determines projection type for grid structure fill
!   \item[latlon\_to\_ij](\ref{latlon_to_ij}) \newline
!     Converts the input lat/lon values to the Cartesian 
!      (i,j) value for the projection type.
!  \end{description}
!EOP

!- Using Map Utils routines to obtain the grid structure 
!   and the grid xpts and ypts information
  type(proj_info) :: proj
  integer :: i
  real, parameter:: RERTH=6.3712E3
  real, parameter:: PI=3.14159265358979,DPR=180./PI
  
  integer :: im, jm,n
  real :: rlat1,rlon1,dx,dy,orient,h
  real :: xmin,xmax,ymin,ymax
  real :: xr,yr

  if(griddesc(1).eq.8) then
     if(gridDesc(13).ne.1 ) then 
        call map_set(PROJ_HRAP,gridDesc(4),gridDesc(5),&
             gridDesc(8)*1000,gridDesc(11),gridDesc(10),0.0,&
             nint(gridDesc(2)),nint(gridDesc(3)),proj)
        do i=1,npts
           call latlon_to_ij(proj,rlat(i),rlon(i),xpts(i),ypts(i))
        enddo
     elseif(gridDesc(13).eq.1) then 
        im=gridDesc(2)
        jm=gridDesc(3)
        rlat1=gridDesc(4)
        rlon1=gridDesc(5)
        orient=gridDesc(7)
        dx=gridDesc(8)
        dy=gridDesc(9)
        if(gridDesc(10).eq.0) then 
           h = 1.0
        else
           h = -1.0
        endif
        
        xmin=0
        xmax=im+1
        ymin=0
        ymax=jm+1
        nret=0
        
        do n=1,npts
           if(abs(rlon(n)).le.360.and.abs(rlat(n)).le.90.and.&
                h*rlat(n).ne.-90) then
              call latlonTohrap(rlat(n),rlon(n),dx,orient,xr,yr)
              
              xpts(n) = xr+(xmax-1)/2+1
              ypts(n) = yr+(ymax-1)/2+1
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
     endif
  endif
end subroutine compute_grid_coord_hrap



