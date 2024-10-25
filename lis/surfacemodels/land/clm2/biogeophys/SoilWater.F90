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
 subroutine SoilWater     (clm, vol_liq, dwat, hk, dhkdw)

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
! Soil hydrology
!
! Method:
! Soil moisture is predicted from a 10-layer model (as with soil 
! temperature), in which the vertical soil moisture transport is governed
! by infiltration, runoff, gradient diffusion, gravity, and root 
! extraction through canopy transpiration.  The net water applied to the
! surface layer is the snowmelt plus precipitation plus the throughfall 
! of canopy dew minus surface runoff and evaporation. 
!
! The vertical water flow in an unsaturated porous media is described by
! Darcy's law, and the hydraulic conductivity and the soil negative 
! potential vary with soil water content and soil texture based on the work 
! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
! integrated over the layer thickness, in which the time rate of change in
! water mass must equal the net flow across the bounding interface, plus the
! rate of internal source or sink. The terms of water flow across the layer
! interfaces are linearly expanded by using first-order Taylor expansion.  
! The equations result in a tridiagonal system equation. 
!
! Note: length units here are all millimeter 
! (in temperature subroutine uses same soil layer 
! structure required but lengths are m)
!
! Richards equation:
!
! d wat      d     d wat d psi
! ----- = - -- [ k(----- ----- - 1) ] + S
!   dt      dz       dz  d wat
!
! where: wat = volume of water per volume of soil (mm**3/mm**3)
! psi = soil matrix potential (mm)
! dt  = time step (s)
! z   = depth (mm)
! dz  = thickness (mm)
! qin = inflow at top (mm h2o /s) 
! qout= outflow at bottom (mm h2o /s)
! s   = source/sink flux (mm h2o /s) 
! k   = hydraulic conductivity (mm h2o /s)
!
!                       d qin                  d qin
! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
!                       d wat(j-1)             d wat(j)
!                ==================|================= 
!                                  < qin 
!
!                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j) 
!
!                                  > qout
!                ==================|================= 
!                        d qout               d qout
! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
!                        d wat(j)             d wat(j+1)
!
!
! Solution: linearize k and psi about d wat and use tridiagonal 
! system of equations to solve for d wat, 
! where for layer j
!
!
! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: SoilWater.F90,v 1.6 2004/11/24 22:56:42 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varpar   , only : nlevsoi
  use clm2_shr_const_mod, only : SHR_CONST_TKFRZ,SHR_CONST_LATICE,SHR_CONST_G
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm !CLM 1-D Module

  real(r8), intent(in) :: vol_liq(1 : nlevsoi) ! soil water per unit volume [mm/mm]

  real(r8), intent(out) :: dwat (1 : nlevsoi)  ! change of soil water [m^3 m-3]
  real(r8), intent(out) :: hk   (1 : nlevsoi)  ! hydraulic conductivity [mm h2o/s]
  real(r8), intent(out) :: dhkdw(1 : nlevsoi)  ! d(hk)/d(vol_liq)

!----Local Variables----------------------------------------------------

  integer j                  ! do loop indices 
  real(r8) amx(1:nlevsoi)    ! "a" left off diagonal of tridiagonal matrix
  real(r8) bmx(1:nlevsoi)    ! "b" diagonal column for tridiagonal matrix
  real(r8) cmx(1:nlevsoi)    ! "c" right off diagonal tridiagonal matrix
  real(r8) z (1 : nlevsoi)   ! layer depth [mm]
  real(r8) dz(1 : nlevsoi)   ! layer thickness [mm]
  real(r8) den               ! used in calculating qin, qout
  real(r8) dqidw0            ! d(qin)/d(vol_liq(i-1))
  real(r8) dqidw1            ! d(qin)/d(vol_liq(i))
  real(r8) dqodw1            ! d(qout)/d(vol_liq(i))
  real(r8) dqodw2            ! d(qout)/d(vol_liq(i+1))
  real(r8) dsmpdw(1:nlevsoi) ! d(smp)/d(vol_liq)
  real(r8) num               ! used in calculating qin, qout
  real(r8) qin               ! flux of water into soil layer [mm h2o/s]
  real(r8) qout              ! flux of water out of soil layer [mm h2o/s]
  real(r8) rmx(1:nlevsoi)    ! "r" forcing term of tridiagonal matrix
  real(r8) s_node            ! soil wetness
  real(r8) s1                ! "s" at interface of layer
  real(r8) s2                ! k*s**(2b+2)
  real(r8) smp(1:nlevsoi)    ! soil matrix potential [mm]
  real(r8) sdamp             ! extrapolates soiwat dependence of evaporation
  real(r8) wimp
!----End Variable List--------------------------------------------------

  wimp = 0.05
  sdamp = 0.
  do j = 1, nlevsoi
     z(j) = clm%z(j)*1.e3
     dz(j) = clm%dz(j)*1.e3
  enddo

!
! Set hydraulic conductivity to zero if effective porosity 5% in any of 
! two neighbor layers or liquid content (theta) less than 0.001
!
  hk = 0.0
  do j = 1, nlevsoi
     if (      (clm%eff_porosity(j) < wimp) &
          .or. (clm%eff_porosity(min(nlevsoi,j+1)) < wimp) &
          .or. (vol_liq(j) <= 1.e-3))then
        hk(j) = 0.
        dhkdw(j) = 0.
     else
        s1 = 0.5*(vol_liq(j)+vol_liq(min(nlevsoi,j+1))) / &
            (0.5*(clm%watsat(j)+clm%watsat(min(nlevsoi,j+1))))
        s2 = clm%hksat(j)*s1**(2.*clm%bsw(j)+2.)
        hk(j) = s1*s2  
        if(hk(j).lt.1e-20) hk(j) = 0.0
        dhkdw(j) = (2.*clm%bsw(j)+3.)*s2*0.5/clm%watsat(j)
        if(j == nlevsoi) dhkdw(j) = dhkdw(j) * 2.
     endif
  enddo
!
! Evaluate hydraulic conductivity, soil matric potential,
! d(smp)/d(vol_liq), and d(hk)/d(vol_liq).
!

  do j = 1, nlevsoi

     if (clm%t_soisno(j)>SHR_CONST_TKFRZ) then

        s_node = max(vol_liq(j)/clm%watsat(j),0.01_r4)
        s_node = min(1.0_r4,s_node)
        smp(j) = -clm%sucsat(j)*s_node**(-clm%bsw(j))
        smp(j) = max(-1.e8, smp(j))        ! Limit soil suction
        dsmpdw(j) = -clm%bsw(j)*smp(j)/(s_node*clm%watsat(j))

     else

!
! When ice is present, the matric potential is only related to temperature
! by (Fuchs et al., 1978: Soil Sci. Soc. Amer. J. 42(3):379-385)
! Unit 1 Joule = 1 (kg m2/s2), (J kg-1)/(m s-2) ==> m ==> 1e3 mm 
!

        smp(j) = 1.e3 * SHR_CONST_LATICE/SHR_CONST_G*(clm%t_soisno(j)-SHR_CONST_TKFRZ)/clm%t_soisno(j)
        smp(j) = max(-1.e8, smp(j))        ! Limit soil suction
        dsmpdw(j) = 0.

     endif
  enddo

!
! Set up r, a, b, and c vectors for tridiagonal solution
!
! Node j=1
!

  j      = 1
  qin    = clm%qflx_infl
  den    = (z(j+1)-z(j))
  num    = (smp(j+1)-smp(j)) - den
  qout   = -hk(j)*num/den
  dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den
  dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den
  rmx(j) =  qin - qout - clm%qflx_tran_veg*clm%rootr(j)
  amx(j) =  0.
  bmx(j) =  dz(j)*(sdamp+1./clm%dtime) + dqodw1
  cmx(j) =  dqodw2

!
! Nodes j=2 to j=nlevsoi-1
!

  do j = 2, nlevsoi - 1
     den    = (z(j) - z(j-1))
     num    = (smp(j)-smp(j-1)) - den
     qin    = -hk(j-1)*num/den
     dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
     dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den
     den    = (z(j+1)-z(j))
     num    = (smp(j+1)-smp(j)) - den
     qout   = -hk(j)*num/den
     dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den
     dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den
     rmx(j) =  qin - qout - clm%qflx_tran_veg*clm%rootr(j)
     amx(j) = -dqidw0
     bmx(j) =  dz(j)/clm%dtime - dqidw1 + dqodw1
     cmx(j) =  dqodw2
  enddo

!
! Node j=nlevsoi
!

  j      = nlevsoi
  den    = (z(j) - z(j-1))
  num    = (smp(j)-smp(j-1)) - den
  qin    = -hk(j-1)*num/den
  dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
  dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den
  qout   =  hk(j)
  dqodw1 =  dhkdw(j)
  rmx(j) =  qin - qout - clm%qflx_tran_veg*clm%rootr(j)
  amx(j) = -dqidw0
  bmx(j) =  dz(j)/clm%dtime - dqidw1 + dqodw1
  cmx(j) =  0.

!
! Solve for dwat
!

  call Tridiagonal (nlevsoi, amx, bmx, cmx, rmx, &
                    dwat)

! 
! Renew the mass of liquid water
!

  do j= 1,nlevsoi
     clm%h2osoi_liq(j) = clm%h2osoi_liq(j) + dwat(j)*dz(j)
  enddo

end subroutine SoilWater
