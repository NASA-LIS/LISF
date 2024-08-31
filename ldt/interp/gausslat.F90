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
! !ROUTINE : gausslat
! \label{gausslat}
!
! !REVISION HISTORY:
!   04-16-92 Mark Iredell; Initial Specification
!   10-20-97 Mark Iredell; Increased precision
!   05-14-02 Urzula Jambor; Reduced limit of eps from e-12 to e-7
!
! !INTERFACE:
subroutine gausslat(jmax,slat,wlat)
  implicit none
! !ARGUMENTS:
  integer       :: jmax
  real          :: slat(jmax)
  real          :: wlat(jmax)
! !DESCRIPTION:
!   This subroutine computes gaussian latitudes
!   Computes cosines of colatitude and gaussian weights
!   on the gaussian latitudes.  the gaussian latitudes are at
!   the zeroes of the legendre polynomial of the given order.
!
!  The arguments are:
!  \begin{description}
!    \item[jmax]
!     input number of latitudes
!    \item[slat]
!     cosines of colatitude
!    \item[wlat]
!     gaussian weights
!  \end{description}
!EOP
  real, parameter :: pi=3.14159265358979
  real, parameter :: eps=1.e-7
  integer, parameter :: jz=50
  real :: c
  integer:: jh, jhe, n, j
  real :: spmax, sp, r
  real :: pk(jmax/2),pkm1(jmax/2),pkm2(jmax/2)
  real :: bz(jz)
  data bz        / 2.4048255577,  5.5200781103, &
       8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679, &
       21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684, &
       33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132, &
       46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550, &
       58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299, &
       71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711, &
       84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819, &
       96.6052679510, 99.7468198587, 102.888374254, 106.029930916, &
       109.171489649, 112.313050280, 115.454612653, 118.596176630, &
       121.737742088, 124.879308913, 128.020877005, 131.162446275, &
       134.304016638, 137.445588020, 140.587160352, 143.728733573, &
       146.870307625, 150.011882457, 153.153458019, 156.295034268 /

  c=(1.-(2./pi)**2)*0.25
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  jh=jmax/2
  jhe=(jmax+1)/2
  r=1./sqrt((jmax+0.5)**2+c)
  do j=1,min(jh,jz)
     slat(j)=cos(bz(j)*r)
  enddo
  do j=jz+1,jh
     slat(j)=cos((bz(jz)+(j-jz)*pi)*r)
  enddo
  spmax=1.
  do while(spmax.gt.eps)
     spmax=0.
     do j=1,jh
        pkm1(j)=1.
        pk(j)=slat(j)
     enddo
     do n=2,jmax
        do j=1,jh
           pkm2(j)=pkm1(j)
           pkm1(j)=pk(j)
           pk(j)=((2*n-1)*slat(j)*pkm1(j)-(n-1)*pkm2(j))/n
        enddo
     enddo
     do j=1,jh
        sp=pk(j)*(1.-slat(j)**2)/(jmax*(pkm1(j)-slat(j)*pk(j)))
        slat(j)=slat(j)-sp
        spmax=max(spmax,abs(sp))
     enddo
  enddo
  do j=1,jh
     wlat(j)=(2.*(1.-slat(j)**2))/(jmax*pkm1(j))**2
     slat(jmax+1-j)=-slat(j)
     wlat(jmax+1-j)=wlat(j)
  enddo
  if(jhe.gt.jh) then
     slat(jhe)=0.
     wlat(jhe)=2./jmax**2
     do n=2,jmax,2
        wlat(jhe)=wlat(jhe)*n**2/(n-1)**2
     enddo
  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  return
end subroutine gausslat
