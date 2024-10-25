!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! Purpose:
!
!    Module to calculate the Gaussian Weights. Public interface is
!    the subroutine "gauaw( a, w, k )".
!
! Method: 
!
!	The algorithm is described in Davis and Rabinowitz,
!      Journal of Research of the NBS, V 56, Jan 1956.
!
! Author: David Williamson, Jim Hack
!
!-----------------------------------------------------------------------
#include "LIS_misc.h"
module gauaw_mod
   save
!
! Public variables
!

   integer, public, parameter :: r16 = selected_real_kind(12)

!
! Public subroutines
!

   public gauaw

!
! Variables private to routines inside this module
!

   real(r16), private ::    pi           ! value of pi
!  real(r16), private, parameter :: one = 1.0     ! 1. in real(r16).  Needed by atan

!
! Functions private to routines inside this module
!

   private bsslzr

contains

   subroutine gauaw(a, w, k)
!-----------------------------------------------------------------------
!
! Calculate sine of latitudes a(k) and weights w(k) for the gaussian
! quadrature. The algorithm is described in Davis and Rabinowitz,
! Journal of Research of the NBS, V 56, Jan 1956.
! The zeros of the bessel function j0, which are obtained from bsslzr,
! are used as a first guess for the abscissa.
!
! Modified 1/23/97 by Jim Rosinski to use real*16 arithmetic in order to
! achieve (nearly) identical weights and latitudes on all machines.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      L. Bath, Jun 1992
!                    L. Buja, Feb 1996
! Reviewed:          D. Williamson, J. Hack, Aug 1992
!                    D. Williamson, J. Hack, Feb 1996
!
!-----------------------------------------------------------------------
!
! $Id: gauaw_mod.F90,v 1.5 2004/05/07 22:18:33 jim Exp $
! $Author: jim $
!
!-----------------------------------------------------------------------
   use LIS_precisionMod
   use clm2_shr_const_mod, only : SHR_CONST_PI
   implicit none
!------------------------------Arguments--------------------------------
!
   integer , intent(in)  :: k      ! number of latitudes pole to pole
   real(r4), intent(out) :: a(k)   ! sine of latitudes
   real(r4), intent(out) :: w(k)   ! gaussian weights
!
!---------------------------Local workspace-----------------------------
!
   real(r16) sinlat(k)    ! sine of latitudes
   real(r16) wgt(k)       ! gaussian weights
 
   real(r16) eps          ! convergence criterion
   real(r16) c            ! constant combination
   real(r16) fk           ! real k
   real(r16) xz           ! abscissa estimate
   real(r16) pkm1         ! |
   real(r16) pkm2         ! |-polynomials
   real(r16) pkmrk        ! |
   real(r16) pk           ! |
   real(r16) sp           ! current iteration latitude increment
   real(r16) avsp         ! |sp|
   real(r16) fn           ! real n
#if ( defined PGF90 )
   parameter (eps = 1.D-15)
#else
   parameter (eps = 1.D-27)
#endif
 
   integer kk           ! k/2 (number of latitudes in hemisphere)
   integer is           ! latitude index
   integer iter         ! iteration counter
   integer n,l          ! indices
!
!-----------------------------------------------------------------------
!
   pi  = SHR_CONST_PI
!
! The value eps, used for convergence tests in the iterations,
! can be changed.  Newton iteration is used to find the abscissas.
!
   c = (1.-(2./pi)**2)*0.25
   fk = k
   kk = k/2
   call bsslzr(sinlat,kk)
   do is=1,kk
     xz = cos(sinlat(is)/sqrt((fk+0.5)**2+c))
!
! This is the first approximation to xz
!
     iter = 0
 10  continue
     pkm2 = 1.
     pkm1 = xz
     iter = iter + 1
     if (iter.gt.10) then
!
! Error exit
!
       write(6,*)'GAUAW:Error exit,no convergence in 10 iterations'
       call endrun
     end if
!
! Computation of the legendre polynomial
!
     do n=2,k
       fn = n
       pk = ((2.*fn-1.)*xz*pkm1-(fn-1.)*pkm2)/fn
       pkm2 = pkm1
       pkm1 = pk
     enddo
     pkm1 = pkm2
     pkmrk = (fk*(pkm1-xz*pk))/(1.-xz**2)
     sp = pk/pkmrk
     xz = xz - sp
     avsp = abs(sp)
     if (avsp.gt.eps) go to 10
     sinlat(is) = xz
     wgt(is) = (2.*(1.-xz**2))/(fk*pkm1)**2
   end do
!
   if (k.ne.kk*2) then
!
! For odd k computation of weight at the equator
!
     sinlat(kk+1) = 0.
     pk = 2./fk**2
     do n=2,k,2
       fn = n
       pk = pk*fn**2/(fn-1.)**2
     end do
     wgt(kk+1) = pk
   end if
!
! Complete the sets of abscissas and weights, using the symmetry.
! Also note truncation from real(r16) to real*8
!
   do n=1,kk
     l = k + 1 - n
     a(n) = sinlat(n)
     a(l) = -sinlat(n)
 
     w(n) = wgt(n)
     w(l) = wgt(n)
   end do
 
   return
   end subroutine gauaw


 
!===========================================================================
 
   subroutine bsslzr(bes, n) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
   use LIS_precisionMod
!
! Return n zeros (or if n>50, approximate zeros), of the Bessel function
! j0,in the array bes. The first 50 zeros will be given exactly, and the
! remaining zeros are computed by extrapolation,and therefore not exact.
!
! Modified 1/23/97 by Jim Rosinski to use real*16 arithmetic
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, D. Williamson, August 1992
! Reviewed:          J. Hack, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer , intent(in) :: n 
   real(r16) , intent(inout) :: bes(n) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer :: j, nn 
   real(r16), dimension(50) :: bz 

   save bz 
!-----------------------------------------------
!------------------------------Arguments--------------------------------
!
!
!---------------------------Local workspace-----------------------------
!
!
   data bz/ 2.4048255577, 5.5200781103, 8.6537279129, 11.7915344391, &
      14.9309177086, 18.0710639679, 21.2116366299, 24.3524715308, &
      27.4934791320, 30.6346064684, 33.7758202136, 36.9170983537, &
      40.0584257646, 43.1997917132, 46.3411883717, 49.4826098974, &
      52.6240518411, 55.7655107550, 58.9069839261, 62.0484691902, &
      65.1899648002, 68.3314693299, 71.4729816036, 74.6145006437, &
      77.7560256304, 80.8975558711, 84.0390907769, 87.1806298436, &
      90.3221726372, 93.4637187819, 96.6052679510, 99.7468198587, &
      102.8883742542, 106.0299309165, 109.1714896498, 112.3130502805, &
      115.4546126537, 118.5961766309, 121.7377420880, 124.8793089132, &
      128.0208770059, 131.1624462752, 134.3040166383, 137.4455880203, &
      140.5871603528, 143.7287335737, 146.8703076258, 150.0118824570, &
      153.1534580192, 156.2950342685/  
!
   nn = n 
   if (n > 50) then 
      bes(50) = bz(50) 
      do j = 51, n 
         bes(j) = bes(j-1) + pi 
      end do 
      nn = 49 
   endif 
   bes(:nn) = bz(:nn) 
   return  
   end subroutine bsslzr 

end module gauaw_mod
