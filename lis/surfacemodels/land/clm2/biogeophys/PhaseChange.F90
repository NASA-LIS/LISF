!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

subroutine PhaseChange (fact,   brr, hs, dhsdT, &   
                        tssbef, xmf, clm ) 

!-----------------------------------------------------------------------
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!-----------------------------------------------------------------------
! Purpose:
! Calculation of the phase change within snow and soil layers:
!
! Method:
! (1) Check the conditions for which the phase change may take place, 
!     i.e., the layer temperature is great than the freezing point 
!     and the ice mass is not equal to zero (i.e. melting), 
!     or the layer temperature is less than the freezing point 
!     and the liquid water mass is not equal to zero (i.e. freezing).
! (2) Assess the rate of phase change from the energy excess (or deficit) 
!     after setting the layer temperature to freezing point.
! (3) Re-adjust the ice and liquid mass, and the layer temperature
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: PhaseChange.F90,v 1.6 2004/11/24 22:56:32 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, only : tfrz, hfus
  use clm2_varpar, only : nlevsoi, nlevsno
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm   ! CLM 1-D Module

  real(r8), intent(in) :: tssbef(clm%snl+1:nlevsoi) ! temperature at previous time step [K]
  real(r8), intent(in) :: brr   (clm%snl+1:nlevsoi) ! 
  real(r8), intent(in) :: fact  (clm%snl+1:nlevsoi) ! temporary variables
  real(r8), intent(in) :: hs                        ! net ground heat flux into the surface
  real(r8), intent(in) :: dhsdT                     ! temperature derivative of "hs"

  real(r8), intent(out) :: xmf                      ! total latent heat of phase change

!----Local Variables----------------------------------------------------

  integer j                        ! do loop index
  real(r8) hm(clm%snl+1:nlevsoi)   ! energy residual [W m-2]
  real(r8) xm(clm%snl+1:nlevsoi)   ! melting or freezing within a time step [kg m-2]
  real(r8) heatr                   ! energy residual or loss after melting or freezing
  real(r8) temp1                   ! temporary variables [kg m-2]
  real(r8), dimension(clm%snl+1:nlevsoi) :: wmass0, wice0!, wliq0
  real(r8)  propor,tinc
!----End Variable List--------------------------------------------------

!
! Initialization 
!
  clm%qflx_snomelt = 0.
  xmf = 0.

  do j = clm%snl+1, nlevsoi
     clm%imelt(j) = 0
     hm(j) = 0.
     xm(j) = 0.
     wice0(j) = clm%h2osoi_ice(j)
!     wliq0(j) = clm%h2osoi_liq(j)
     wmass0(j) = clm%h2osoi_ice(j) + clm%h2osoi_liq(j)
  enddo

!
! Melting identification
! If ice exists above melt point, melt some to liquid.
!
  do j = clm%snl+1, nlevsoi
     if (clm%h2osoi_ice(j) > 0. .AND. clm%t_soisno(j) > tfrz) then
        clm%imelt(j) = 1
        clm%t_soisno(j) = tfrz
     endif

!
! Freezing identification
! If liquid exists below melt point, freeze some to ice.
!

     if (clm%h2osoi_liq(j) > 0. .AND. clm%t_soisno(j) < tfrz) then
        clm%imelt(j) = 2
        clm%t_soisno(j) = tfrz
     endif
  enddo

!
! If snow exists, but its thickness is less than the critical value (0.01 m)
!

  if (clm%snl+1 == 1 .AND. clm%h2osno > 0.) then
     if (clm%t_soisno(1) > tfrz) then
        clm%imelt(1) = 1
        clm%t_soisno(1) = tfrz
     endif
  endif

!
! Calculate the energy surplus and loss for melting and freezing
!

  do j = clm%snl+1, nlevsoi
     if (clm%imelt(j) > 0) then
        tinc = clm%t_soisno(j)-tssbef(j)
        if (j > clm%snl+1) then
           hm(j) = brr(j) - tinc/fact(j) 
        else
           hm(j) = hs + dhsdT*tinc + brr(j) - tinc/fact(j) 
        endif
     endif
  enddo

!
! These two errors were checked carefully.  They result from the 
! computed error of "Tridiagonal-Matrix" in subroutine "thermal".
!

  do j = clm%snl+1, nlevsoi
     if (clm%imelt(j) == 1 .AND. hm(j) < 0.) then
        hm(j) = 0.
        clm%imelt(j) = 0
     endif

     if (clm%imelt(j) == 2 .AND. hm(j) > 0.) then
        hm(j) = 0.
        clm%imelt(j) = 0
     endif
  enddo

!
! The rate of melting and freezing
!

  do j = clm%snl+1, nlevsoi

     if (clm%imelt(j) > 0 .and. abs(hm(j)) > .0) then
        xm(j) = hm(j)*clm%dtime/hfus                        ! kg m-2

!
! If snow exists, but its thickness is less than the critical value
! (1 cm). Note: more work is needed to determine how to tune the
! snow depth for this case
!

        if (j == 1) then
           if ((clm%snl+1 == 1) .AND. (clm%h2osno > 0.) .AND. (xm(j) > 0.)) then
              temp1 = clm%h2osno                                        ! kg m-2
              clm%h2osno = max(0._r4,temp1-xm(j))
              propor = clm%h2osno/temp1
              clm%snowdp = propor * clm%snowdp
              heatr = hm(j) - hfus*(temp1-clm%h2osno)/clm%dtime         ! W m-2
              if (heatr > 0.) then
                 xm(j) = heatr*clm%dtime/hfus                           ! kg m-2
                 hm(j) = heatr                                          ! W m-2
              else
                 xm(j) = 0.
                 hm(j) = 0.
              endif
              clm%qflx_snomelt = max(0._r4,(temp1-clm%h2osno))/clm%dtime   ! kg/(m2 s)
              xmf = hfus*clm%qflx_snomelt
           endif
        endif

        heatr = 0.
        if (xm(j) > 0.) then
           clm%h2osoi_ice(j) = max(0._r4, wice0(j)-xm(j))
           heatr = hm(j) - hfus*(wice0(j)-clm%h2osoi_ice(j))/clm%dtime
        else if (xm(j) < 0.) then
           clm%h2osoi_ice(j) = min(wmass0(j), wice0(j)-xm(j))
           heatr = hm(j) - hfus*(wice0(j)-clm%h2osoi_ice(j))/clm%dtime  
        endif

        clm%h2osoi_liq(j) = max(0._r4,wmass0(j)-clm%h2osoi_ice(j))

        if (abs(heatr) > 0.) then
           if (j > clm%snl+1) then
              clm%t_soisno(j) = clm%t_soisno(j) + fact(j)*heatr
           else
              clm%t_soisno(j) = clm%t_soisno(j) + fact(j)*heatr/(1.-fact(j)*dhsdT)
           endif
           if (clm%h2osoi_liq(j)*clm%h2osoi_ice(j)>0.) clm%t_soisno(j) = tfrz
        endif

        xmf = xmf + hfus * (wice0(j)-clm%h2osoi_ice(j))/clm%dtime

        if (clm%imelt(j) == 1 .AND. j < 1) then
           clm%qflx_snomelt = clm%qflx_snomelt + max(0._r4,(wice0(j)-clm%h2osoi_ice(j)))/clm%dtime  
        endif

     endif

  enddo

!
! needed for output to history file
!

!  clm%eflx_snomelt = clm%qflx_snomelt * hfus  

end subroutine PhaseChange
