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
! !ROUTINE: compute_grid_coord_gauss
! \label{compute_grid_coord_gauss}
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
subroutine compute_grid_coord_gauss(gridDesc,npts,fill,xpts,ypts,&
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
!  the specified domain for a gaussian cylindrical projection.
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
!   \item[gausslat](\ref{gausslat} \newline
!     Computes latitude values in gaussian
!  \end{description}
!EOP
  real, parameter :: pi=3.14159265358979
  integer, parameter :: jgmax=2000
  real :: dpr
  integer :: im,jm, jg, ja, n
  real :: rlat1,rlon1, rlat2, rlon2
  real :: hi, wb
  real :: dlon
  real :: xmin,xmax,ymin,ymax
  real :: alat(0:jgmax+1),blat(jgmax)
  integer :: iscan,jscan,nscan, iret
  real :: yptsa, yptsb
  integer :: jh, j1, j2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  dpr =180./pi 
  if(gridDesc(1).eq.4.and.gridDesc(10)*2.le.jgmax) then
     im=gridDesc(2)
     jm=gridDesc(3)
     rlat1=gridDesc(4)
     rlon1=gridDesc(5)
     rlat2=gridDesc(7)
     rlon2=gridDesc(8)
     jg=gridDesc(10)*2
     iscan=mod(nint(gridDesc(11))/128,2)
     jscan=mod(nint(gridDesc(11))/64,2)
     nscan=mod(nint(gridDesc(11))/32,2)
     hi=(-1.)**iscan
     jh=(-1)**jscan
     dlon=hi*(mod(hi*(rlon2-rlon1)-1+3600,360.)+1)/(im-1)
     call gausslat(jg,alat(1),blat)
     do ja=1,jg
        alat(ja)=dpr*asin(alat(ja))
     enddo
     alat(0)=180.-alat(1)
     alat(jg+1)=-alat(0)
     j1=1
     do while(j1.lt.jg.and.rlat1.lt.(alat(j1)+alat(j1+1))/2)
        j1=j1+1
     enddo
     j2=j1+jh*(jm-1)
     xmin=0
     xmax=im+1
     if(im.eq.nint(360/abs(dlon))) xmax=im+2
     ymin=0.5
     ymax=jm+0.5
     nret=0
! translate grid coordinates to earth coordinates

     if(abs(dlon-gridDesc(9)).gt.0.01) then
        print*, 'problem with the domain calculations : gdswiz04'
        stop
     endif
     do n=1,npts
        xpts(n)=fill
        ypts(n)=fill
        if(abs(rlon(n)).le.360.and.abs(rlat(n)).le.90) then
           xpts(n)=1+hi*mod(hi*(rlon(n)-rlon1)+3600,360.)/dlon
           ja=min(int((jg+1)/180.*(90-rlat(n))),jg)
           if(rlat(n).gt.alat(ja)) ja=max(ja-2,0)
           if(rlat(n).lt.alat(ja+1)) ja=min(ja+2,jg)
           if(rlat(n).gt.alat(ja)) ja=ja-1
           if(rlat(n).lt.alat(ja+1)) ja=ja+1
           yptsa=1+jh*(ja-j1)
           yptsb=1+jh*(ja+1-j1)
           wb=(alat(ja)-rlat(n))/(alat(ja)-alat(ja+1))
           ypts(n)=yptsa+wb*(yptsb-yptsa)
           if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. & 
                ypts(n).ge.ymin.and.ypts(n).le.ymax) then
              nret=nret+1
           else
              xpts(n)=fill
              ypts(n)=fill
           endif
        endif
     enddo

! projection unrecognized
  else
     iret=-1
     do n=1,npts
        xpts(n)=fill
        ypts(n)=fill
     enddo     
  endif
end subroutine compute_grid_coord_gauss
