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
! !ROUTINE: compute_earth_coord_gauss
! \label{compute_earth_coord_gauss}
!
! !REVISION HISTORY: 
!   04-10-96 Mark Iredell;  Initial Specification
!   05-27-04 Sujay Kumar; Modified verision with floating point arithmetic. 
!
! !INTERFACE:
subroutine compute_earth_coord_gauss(gridDesc,npts,fill,xpts,ypts,&
     rlon,rlat,nret)

  implicit none
! !ARGUMENTS: 
  real            :: gridDesc(20)
  integer         :: npts
  real            :: fill
  real            :: xpts(npts),ypts(npts)
  real            :: rlat(npts)
  real            :: rlon(npts)

! !DESCRIPTION:
!  This subroutine computes the earth coordinates of 
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
!     grid x point coordinates
!    \item[ypts]
!     grid y point coordinates
!    \item[rlat]    
!     output latitudes in degrees
!    \item[rlon]    
!     output longitudes in degrees
!    \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[gausslat](\ref{gausslat}) \newline
!     Computes latitude values in gaussian
!  \end{description}
!EOP

  real, parameter :: pi=3.14159265358979
  integer, parameter :: jgmax=2000
  real :: dpr
  real ::  rlata,rlatb  
  integer :: nret
  integer :: im,jm, jg, j, ja, n
  real :: rlat1,rlon1, rlat2, rlon2
  real :: hi, wb
  real :: dlon
  real :: xmin,xmax,ymin,ymax
  real :: alat(0:jgmax+1),blat(jgmax)
  integer :: iscan,jscan,nscan, iret
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
     do n=1,npts
        if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. & 
             ypts(n).ge.ymin.and.ypts(n).le.ymax) then
           rlon(n)=mod(rlon1+dlon*(xpts(n)-1)+3600,360.)
           j=min(int(ypts(n)),jm)
           rlata=alat(j1+jh*(j-1))
           rlatb=alat(j1+jh*j)
           wb=ypts(n)-j
           rlat(n)=rlata+wb*(rlatb-rlata)
           nret=nret+1
        else
           rlon(n)=fill
           rlat(n)=fill
        endif
     enddo

! projection unrecognized
  else
     iret=-1
     
     do n=1,npts
        rlon(n)=fill
        rlat(n)=fill
     enddo
  endif
end subroutine compute_earth_coord_gauss
