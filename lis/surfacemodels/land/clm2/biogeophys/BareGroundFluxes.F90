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

subroutine BareGroundFluxes (tg,     thm,   qg,    thv,   z0mg,   &
                            z0hg,   z0qg,  dqgdT, htvp,  beta,   &
                            zii,    ur,    dlrad, ulrad, cgrnds, &
                            cgrndl, cgrnd, clm    )

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
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! Method:
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: BareGroundFluxes.F90,v 1.6 2004/11/24 22:56:16 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, only : cpair, vkc, grav
  use clm2_shr_const_mod, only : SHR_CONST_RGAS
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm !CLM 1-D Module

  real(r8), intent(in) :: tg    ! ground surface temperature [K]
  real(r8), intent(in) :: thm   ! intermediate variable (forc_t+0.0098*forc_hgt_t)
  real(r8), intent(in) :: qg    ! specific humidity at ground surface [kg kg-1]
  real(r8), intent(in) :: thv   ! virtual potential temperature (kelvin)
  real(r8), intent(in) :: z0mg  ! roughness length, momentum [m]
  real(r8), intent(in) :: dqgdT ! temperature derivative of "qg"
  real(r8), intent(in) :: htvp  ! latent heat of evaporation (/sublimation) [J kg-1]
  real(r8), intent(in) :: beta  ! coefficient of conective velocity [-]
  real(r8), intent(in) :: zii   ! convective boundary height [m]
  real(r8), intent(in) :: ur    ! wind speed at reference height [m s-1]

  real(r8), intent(inout) :: z0hg   ! roughness length, sensible heat [m]
  real(r8), intent(inout) :: z0qg   ! roughness length, latent heat [m]
  real(r8), intent(inout) :: cgrnd  ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
  real(r8), intent(inout) :: cgrndl ! deriv. of soil sensible heat flux wrt soil temp
                                    ! [w/m2/k]
  real(r8), intent(inout) :: cgrnds ! deriv. of soil latent heat flux wrt soil temp
                                    ! [w/m**2/k]

  real(r8), intent(out) :: dlrad ! downward longwave radiation below the canopy [W m-2]
  real(r8), intent(out) :: ulrad ! upward longwave radiation above the canopy [W m-2]

!----Local Variables----------------------------------------------------

  integer nmozsgn  ! number of times moz changes sign
  integer niters   ! maximum number of iterations for surface temperature
  integer iter     ! iteration index
  real(r8) zldis   ! reference height "minus" zero displacement height [m]
  real(r8) displa  ! displacement height [m]
  real(r8) zeta    ! dimensionless height used in Monin-Obukhov theory
  real(r8) wc      ! convective velocity [m s-1]
  real(r8) dth     ! diff of virtual temp. between ref. height and surface 
  real(r8) dthv    ! diff of vir. poten. temp. between ref. height and surface
  real(r8) dqh     ! diff of humidity between ref. height and surface
  real(r8) obu     ! Monin-Obukhov length (m)
  real(r8) um      ! wind speed including the stablity effect [m s-1]
  real(r8) temp1   ! relation for potential temperature profile
  real(r8) temp2   ! relation for specific humidity profile
  real(r8) temp3   ! relation for 2m potentail temperture profile
  real(r8) temp4   ! relation for 2m specific humidity profile
  real(r8) ustar   ! friction velocity [m s-1]
  real(r8) tstar   ! temperature scaling parameter
  real(r8) qstar   ! moisture scaling parameter
  real(r8) thvstar ! virtual potential temperature scaling parameter
  real(r8) cf      ! heat transfer coefficient from leaves [-]
  real(r8) ram     ! aerodynamical resistance [s/m]
  real(r8) rah     ! thermal resistance [s/m]
  real(r8) raw     ! moisture resistance [s/m]
  real(r8) raih    ! temporary variable [kg m-2/s]
  real(r8) raiw    ! temporary variable [kg m-2/s]
  real(r8) obuold  ! monin-obukhov length from previous iteration

!----End Variable List--------------------------------------------------

!
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!

!
! Initialization variables
!

     dlrad  = 0.
     ulrad  = 0.

     nmozsgn = 0
     obuold = 0.
     dth   = thm-tg
     dqh   = clm%forc_q-qg
     dthv  = dth*(1.+0.61*clm%forc_q)+0.61*clm%forc_th*dqh
     zldis = clm%forc_hgt_u-0.

!
! Initialize Monin-Obukhov length and wind speed
!

     call MoninObukIni(ur, thv, dthv, zldis, z0mg, &
                       um, obu  )

!
! Begin stability iteration
! Determine friction velocity, and potential temperature and humidity
! profiles of the surface boundary layer
!

     niters=3
     do iter = 1, niters

        displa = 0.0_r4
        call FrictionVelocity(displa,z0mg,z0hg,z0qg,obu, &
                              iter,ur,um,ustar,temp1,temp2, &
                              temp3,temp4,clm)

        tstar = temp1*dth
        qstar = temp2*dqh
        z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
        z0qg = z0hg

        thvstar=tstar*(1.+0.61*clm%forc_q) + 0.61*clm%forc_th*qstar
        zeta=zldis*vkc*grav*thvstar/(ustar**2*thv)
        if (zeta >= 0.) then     !stable
           zeta = min(2._r4,max(zeta,0.01_r4))
           um = max(ur,0.1_r4)
        else                     !unstable
           zeta = max(-100._r4,min(zeta,-0.01_r4))
           wc = beta*(-grav*ustar*thvstar*zii/thv)**0.333
           um = sqrt(ur*ur+wc*wc)
        endif
        obu = zldis/zeta

        if (obuold*obu < 0.) nmozsgn = nmozsgn+1
        if (nmozsgn >= 4) EXIT

        obuold = obu

     enddo                       ! end stability iteration

!
! Determine aerodynamic resistances
!

     clm%ch = ustar*ustar/um ! Add-in for ALMA output
     clm%chs2 = temp3*ustar ! For calc of 2M T
     clm%cqs2 = temp4*ustar ! For calc of 2M Q

     ram    = 1./(ustar*ustar/um)
     rah    = 1./(temp1*ustar)
     raw    = 1./(temp2*ustar)
     raih   = clm%forc_rho*cpair/rah
     raiw   = clm%forc_rho/raw

!
! Get derivative of fluxes with respect to ground temperature
!

     cgrnds = raih
     cgrndl = raiw*dqgdT
     cgrnd  = cgrnds + htvp*cgrndl

!
! Surface fluxes of momentum, sensible and latent heat
! using ground temperatures from previous time step
!

     clm%taux   = -clm%forc_rho*clm%forc_u/ram
     clm%tauy   = -clm%forc_rho*clm%forc_v/ram
     clm%eflx_sh_grnd  = -raih*dth
     clm%qflx_evap_soi  = -raiw*dqh
     clm%eflx_sh_tot  = clm%eflx_sh_grnd
     clm%qflx_evap_tot  = clm%qflx_evap_soi

!
! 2 m height air temperature
!

     clm%t_ref2m = (tg+temp1*dth * 1./vkc *log((2.+z0hg)/z0hg))

!
! Variables needed by history tape
!
!     print*, 'before..',clm%t_veg
     clm%t_veg = clm%forc_t
!     print*, 'after..',clm%t_veg
!     clm%btran = 0.
     clm%rootr(:) = 0.
     cf = clm%forc_pbot/(SHR_CONST_RGAS*0.001*thm)*1.e06
     clm%rssun = 1./1.e15 * cf
     clm%rssha = 1./1.e15 * cf

end subroutine BareGroundFluxes
