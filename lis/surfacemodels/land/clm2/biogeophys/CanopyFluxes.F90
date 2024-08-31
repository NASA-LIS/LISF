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

subroutine CanopyFluxes (z0mv,   z0hv,  z0qv,  thm,   th,     &
                         thv,    tg,    qg,    dqgdT, htvp,   &
                         emv,    emg,   dlrad, ulrad, cgrnds, &
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
! This subroutine:
! 1. Calculates the leaf temperature: 
! 2. Calculates the leaf fluxes, transpiration, photosynthesis and 
!    updates the dew accumulation due to evaporation.
!
! Method:
! Use the Newton-Raphson iteration to solve for the foliage 
! temperature that balances the surface energy budget:
!
! f(t_veg) = Net radiation - Sensible - Latent = 0
! f(t_veg) + d(f)/d(t_veg) * dt_veg = 0     (*)
!
! Note:
! (1) In solving for t_veg, t_grnd is given from the previous timestep.
! (2) The partial derivatives of aerodynamical resistances, which cannot 
!     be determined analytically, are ignored for d(H)/dT and d(LE)/dT
! (3) The weighted stomatal resistance of sunlit and shaded foliage is used 
! (4) Canopy air temperature and humidity are derived from => Hc + Hg = Ha
!                                                          => Ec + Eg = Ea
! (5) Energy loss is due to: numerical truncation of energy budget equation
!     (*); and "ecidif" (see the code) which is dropped into the sensible 
!     heat 
! (6) The convergence criteria: the difference, del = t_veg(n+1)-t_veg(n) and 
!     del2 = t_veg(n)-t_veg(n-1) less than 0.01 K, and the difference of 
!     water flux from the leaf between the iteration step (n+1) and (n) 
!     less than 0.1 W m-2; or the iterative steps over 40.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: CanopyFluxes.F90,v 1.8 2004/11/24 22:56:20 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, only : sb, cpair, hvap, vkc, grav, denice, denh2o, tfrz
  use clm2_varpar, only : nlevsoi

  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm !CLM 1-D Module

  real(r8), intent(in) :: z0mv  ! roughness length, momentum [m]
  real(r8), intent(in) :: z0hv  ! roughness length, sensible heat [m]
  real(r8), intent(in) :: z0qv  ! roughness length, latent heat [m]
  real(r8), intent(in) :: thm   ! intermediate variable (forc_t+0.0098*forc_hgt_t)
  real(r8), intent(in) :: th    ! potential temperature (kelvin)
  real(r8), intent(in) :: thv   ! virtual potential temperature (kelvin)
  real(r8), intent(in) :: tg    ! ground surface temperature [K]
  real(r8), intent(in) :: qg    ! specific humidity at ground surface [kg kg-1]
  real(r8), intent(in) :: dqgdT ! temperature derivative of "qg"
  real(r8), intent(in) :: htvp  ! latent heat of evaporation (/sublimation) [J kg-1]
  real(r8), intent(in) :: emv   ! ground emissivity
  real(r8), intent(in) :: emg   ! vegetation emissivity

  real(r8), intent(inout) :: cgrnd  ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
  real(r8), intent(inout) :: cgrndl ! deriv, of soil sensible heat flux wrt soil temp 
                                    ! [w/m2/k]
  real(r8), intent(inout) :: cgrnds ! deriv of soil latent heat flux wrt soil temp 
                                    ! [w/m2/k]

  real(r8), intent(out) :: dlrad    ! downward longwave radiation below the canopy [W m-2]
  real(r8), intent(out) :: ulrad    ! upward longwave radiation above the canopy [W m-2]
  real(r8) :: parsha
  real(r8) :: btran
  real(r8) :: csoilc
!----Local Variables----------------------------------------------------

  real(r8) zldis       ! reference height "minus" zero displacement height [m]
  real(r8) zii         ! convective boundary layer height [m]
  real(r8) zeta        ! dimensionless height used in Monin-Obukhov theory
  real(r8) beta        ! coefficient of conective velocity [-]
  real(r8) wc          ! convective velocity [m s-1]
  real(r8) dth         ! diff of virtual temp. between ref. height and surface 
  real(r8) dthv        ! diff of vir. poten. temp. between ref. height and surface
  real(r8) dqh         ! diff of humidity between ref. height and surface
  real(r8) obu         ! Monin-Obukhov length (m)
  real(r8) um          ! wind speed including the stablity effect [m s-1]
  real(r8) ur          ! wind speed at reference height [m s-1]
  real(r8) uaf         ! velocity of air within foliage [m s-1]
  real(r8) temp1       ! relation for potential temperature profile
  real(r8) temp2       ! relation for specific humidity profile
  real(r8) temp3       ! relation for 2m potential temperature profile
  real(r8) temp4       ! relation for 2m specific humidity profile
  real(r8) ustar       ! friction velocity [m s-1]
  real(r8) tstar       ! temperature scaling parameter
  real(r8) qstar       ! moisture scaling parameter
  real(r8) thvstar     ! virtual potential temperature scaling parameter
  real(r8) taf         ! air temperature within canopy space [K]
  real(r8) qaf         ! humidity of canopy air [kg kg-1]
  real(r8) rpp         ! fraction of potential evaporation from leaf [-]
  real(r8) rppdry      ! fraction of potential evaporation through transp [-]
  real(r8) cf          ! heat transfer coefficient from leaves [-]
  real(r8) rb          ! leaf boundary layer resistance [s/m]
  real(r8) ram(2)      ! aerodynamical resistance [s/m]
  real(r8) rah(2)      ! thermal resistance [s/m]
  real(r8) raw(2)      ! moisture resistance [s/m]
  real(r8) wta         ! heat conductance for air [m s-1]
  real(r8) wtg         ! heat conductance for ground [m s-1]
  real(r8) wtl         ! heat conductance for leaf [m s-1]
  real(r8) wta0        ! normalized heat conductance for air [-]
  real(r8) wtl0        ! normalized heat conductance for leaf [-]
  real(r8) wtg0        ! normalized heat conductance for ground [-]
  real(r8) wtal        ! normalized heat conductance for air and leaf [-]
  real(r8) wtgl        ! normalized heat conductance for leaf and ground [-]
  real(r8) wtga        ! normalized heat cond. for air and ground  [-]
  real(r8) wtaq        ! latent heat conductance for air [m s-1]
  real(r8) wtlq        ! latent heat conductance for leaf [m s-1]
  real(r8) wtgq        ! latent heat conductance for ground [m s-1]
  real(r8) wtaq0       ! normalized latent heat conductance for air [-]
  real(r8) wtlq0       ! normalized latent heat conductance for leaf [-]
  real(r8) wtgq0       ! normalized heat conductance for ground [-]
  real(r8) wtalq       ! normalized latent heat cond. for air and leaf [-]
  real(r8) wtglq       ! normalized latent heat cond. for leaf and ground [-]
  real(r8) wtgaq       ! normalized latent heat cond. for air and ground [-]
  real(r8) el          ! vapor pressure on leaf surface [pa]
  real(r8) deldT       ! derivative of "el" on "t_veg" [pa/K]
  real(r8) qsatl       ! leaf specific humidity [kg kg-1]
  real(r8) qsatldT     ! derivative of "qsatl" on "t_veg"
  real(r8) air,bir,cir ! atmos. radiation temporay set
  real(r8) dc1,dc2     ! derivative of energy flux [W m-2/K]
  real(r8) delt        ! temporary
  real(r8) delq        ! temporary
  real(r8) del         ! absolute change in leaf temp in current iteration [K]
  real(r8) del2        ! change in leaf temperature in previous iteration [K]
  real(r8) dele        ! change in latent heat flux from leaf [K]
  real(r8) delmax      ! maximum change in  leaf temperature [K]
  real(r8) dels        ! change in leaf temperature in current iteration [K]
  real(r8) det         ! maximum leaf temp. change in two consecutive iter [K]
  real(r8) dlemin      ! max limit for energy flux convergence [w/m2]
  real(r8) dtmin       ! max limit for temperature convergence [K]
  real(r8) efeb        ! latent heat flux from leaf (previous iter) [mm s-1]
  real(r8) efeold      ! latent heat flux from leaf (previous iter) [mm s-1]
  real(r8) efpot       ! potential latent energy flux [kg m-2/s]
  real(r8) efe         ! water flux from leaf [mm s-1]
  real(r8) efsh        ! sensible heat from leaf [mm s-1]
  real(r8) obuold      ! monin-obukhov length from previous iteration
  real(r8) tlbef       ! leaf temperature from previous iteration [K]
  real(r8) ecidif      ! excess energies [W m-2]
  real(r8) err         ! balance error
  real(r8) erre        ! balance error
  real(r8) co2         ! atmospheric co2 concentration (pa)
  real(r8) o2          ! atmospheric o2 concentration (pa)
  real(r8) svpts       ! saturation vapor pressure at t_veg (pa)
  real(r8) eah         ! canopy air vapor pressure (pa)
  real(r8) s_node      ! vol_liq/eff_porosity
  real(r8) smp_node    ! matrix potential
  real(r8) vol_ice(1:nlevsoi)      ! partial volume of ice lens in layer
  real(r8) eff_porosity(1:nlevsoi) ! effective porosity in layer
  real(r8) vol_liq(1:nlevsoi)      ! partial volume of liquid water in layer
  real(r8) rresis(1:nlevsoi)       ! soil water contribution to root resistance

! Constant atmospheric co2 and o2
  real(r8) po2                   ! partial pressure  o2 (mol/mol)
  real(r8) pco2                  ! partial pressure co2 (mol/mol)
  data po2,pco2 /0.209,355.e-06/

  real(r8) :: mpe = 1.e-6        ! prevents overflow error if division by zero

  integer i       ! loop index
  integer itlef   ! counter for leaf temperature iteration [-]
  integer itmax   ! maximum number of iteration [-]
  integer itmin   ! minimum number of iteration [-]
  integer nmozsgn ! number of times stability changes sign
  real(r8):: dt_veg
  real(r8) :: smpmax
!----End Variable List--------------------------------------------------

!
! Initialization
!
  smpmax = -1.5e5
  csoilc = 0.004

  del   = 0.0  ! change in leaf temperature from previous iteration
  itlef = 0    ! counter for leaf temperature iteration
  efeb  = 0.0  ! latent head flux from leaf for previous iteration

  wtlq = 0.0
  wtlq0 = 0.0
  wtgq0 = 0.0
  wtalq = 0.0
  wtgaq = 0.0
  wtglq = 0.0
  wtaq = 0.0
  wtgq = 0.0
  wtaq0 = 0.0
  wtlq0 = 0.0
  wtgq0 = 0.0
  wtalq = 0.0
  wtgaq = 0.0
  wtglq = 0.0

!
! Assign iteration parameters
!

  delmax = 1.0  ! maximum change in  leaf temperature
  itmax  = 40   ! maximum number of iteration
  itmin  = 2    ! minimum number of iteration
  dtmin  = 0.01 ! max limit for temperature convergence
  dlemin = 0.1  ! max limit for energy flux convergence

!
! Effective porosity of soil, partial volume of ice and liquid (needed
! for btran)
!

  do i = 1,nlevsoi
     vol_ice(i) = min(clm%watsat(i), clm%h2osoi_ice(i)/(clm%dz(i)*denice))
     eff_porosity(i) = clm%watsat(i)-vol_ice(i)
     vol_liq(i) = min(eff_porosity(i), clm%h2osoi_liq(i)/(clm%dz(i)*denh2o))
  enddo

!
! Root resistance factors
!

  btran = 1.e-10
  do i = 1,nlevsoi
     if (clm%t_soisno(i) > tfrz) then
        s_node = max(vol_liq(i)/eff_porosity(i),0.01_r4)
        smp_node = max(smpmax, -clm%sucsat(i)*s_node**(-clm%bsw(i)))
        rresis(i) = (1.-smp_node/smpmax)/(1.+clm%sucsat(i)/smpmax)
        clm%rootr(i) = clm%rootfr(i)*rresis(i)
        btran = btran + clm%rootr(i)
     else
        clm%rootr(i) = 0.
     endif
  enddo

!
! Normalize root resistances to get layer contribution to ET
!

  do i = 1,nlevsoi
     clm%rootr(i)  = clm%rootr(i)/btran
  enddo

!
! Net absorbed longwave radiation by canopy and ground
! =air+bir*t_veg**4+cir*t_grnd**4
!

  air =   emv * (1.+(1.-emv)*(1.-emg)) * clm%forc_lwrad
  bir = - (2.-emv*(1.-emg)) * emv * sb
  cir =   emv*emg*sb

!
! Saturated vapor pressure, specific humidity, and their derivatives
! at the leaf surface
!

  call QSatclm (clm%t_veg, clm%forc_pbot, el, deldT, qsatl, &
             qsatldT)

!
! Determine atmospheric co2 and o2
!

  co2 = pco2*clm%forc_pbot
  o2  = po2*clm%forc_pbot

!
! Initialize flux profile
!

  nmozsgn = 0
  obuold = 0.
  zii=1000.         ! m  (pbl height)
  beta=1.           ! -  (in computing W_*)

  taf = (tg + thm)/2.
  qaf = (clm%forc_q+qg)/2.

  ur = max(1.0_r8,sqrt(clm%forc_u*clm%forc_u+clm%forc_v*clm%forc_v))
  dth = thm-taf
  dqh = clm%forc_q-qaf
  dthv = dth*(1.+0.61*clm%forc_q)+0.61*th*dqh

  zldis = clm%forc_hgt_u - clm%displa

!
! Initialize Monin-Obukhov length and wind speed
!
  call MoninObukIni(ur, thv, dthv, zldis, z0mv, &
                    um, obu  )

!
! Begin stability iteration
!

  ITERATION : do while (itlef <= itmax) 

     tlbef = clm%t_veg
     del2 = del

!
! Determine friction velocity, and potential temperature and humidity
! profiles of the surface boundary layer
!


     call FrictionVelocity (clm%displa, z0mv,  z0hv,  z0qv,  obu, &
                            itlef+1, ur, um, ustar, temp1, temp2, &
                            temp3,temp4,clm)
    
!
! Determine aerodynamic resistances
!

!=== LDAS modification

      clm%ch = ustar*ustar/um ! Add-in for ALMA output
      clm%chs2 = temp3*ustar ! For calc of 2M T
      clm%cqs2 = temp4*ustar ! For calc of 2M Q

     ram(1)=1./(ustar*ustar/um)
     rah(1)=1./(temp1*ustar) 
     raw(1)=1./(temp2*ustar) 
     
!
! Bulk boundary layer resistance of leaves
!

     uaf = um*sqrt( 1./(ram(1)*um) )
     cf = 0.01/(sqrt(uaf)*sqrt(clm%dleaf))
     rb = 1./(cf*uaf)

!
! Aerodynamic resistances raw and rah between heights zpd+z0h and z0hg.
! if no vegetation, rah(2)=0 because zpd+z0h = z0hg.
! (Dickinson et al., 1993, pp.54)
!

     ram(2) = 0.               ! not used
     rah(2) = 1./(csoilc*uaf)
     raw(2) = rah(2) 

!
! Stomatal resistances for sunlit and shaded fractions of canopy.
! Done each iteration to account for differences in eah, tv.
!
     parsha = 0._r4 
     svpts = el                        ! pa
     eah = clm%forc_pbot * qaf / 0.622 ! pa

     call Stomata(mpe,        clm%parsun, svpts,     eah,      thm,        &
                  o2,         co2,        btran, rb,       clm%rssun,  &
                  clm%psnsun, clm%qe25,   clm%vcmx25,clm%mp,   clm%c3psn,  &
                  clm       )
                                                      
     call Stomata(mpe,        parsha, svpts,     eah,      thm,        &
                  o2,         co2,        btran, rb,       clm%rssha,  &
                  clm%psnsha, clm%qe25,   clm%vcmx25,clm%mp,   clm%c3psn,  &
                  clm       )

!
! Heat conductance for air, leaf and ground  
!

     call SensibleHCond(rah(1), rb,   rah(2), wta,  wtl,  &
                        wtg,    wta0, wtl0,   wtg0, wtal, &
                        wtga,   wtgl, clm     )

!
! Fraction of potential evaporation from leaf
!

     if (clm%fdry .gt. 0.0) then
        rppdry  = clm%fdry*rb*(clm%laisun/(rb+clm%rssun) + clm%laisha/(rb+clm%rssha))/ &
                  clm%elai
     else
        rppdry = 0.0
     endif
     efpot = clm%forc_rho*wtl*(qsatl-qaf)

     if (efpot > 0.) then
       if (btran > 1.e-10) then

        clm%qflx_tran_veg = efpot*rppdry
        rpp = rppdry + clm%fwet

!
! No transpiration if btran below 1.e-10
!

       else
        rpp = clm%fwet
        clm%qflx_tran_veg = 0.
       endif

!
! Check total evapotranspiration from leaves
!

       rpp = min(rpp, (clm%qflx_tran_veg+clm%h2ocan/clm%dtime)/efpot)

     else

!
! No transpiration if potential evaporation less than zero
!

       rpp = 1.
       clm%qflx_tran_veg = 0.

     endif

!
! Update conductances for changes in rpp 
! Latent heat conductances for ground and leaf.
! Air has same conductance for both sensible and latent heat.
!

     call LatentHCond(raw(1), rb,    raw(2), rpp,   wtaq,  &
                      wtlq,   wtgq,  wtaq0,  wtlq0, wtgq0, &
                      wtalq,  wtgaq, wtglq,  clm    ) 

     dc1 = clm%forc_rho*cpair*wtl
     dc2 = hvap*clm%forc_rho*wtlq

     efsh = dc1*(wtga*clm%t_veg-wtg0*tg-wta0*thm)
     efe = dc2*(wtgaq*qsatl-wtgq0*qg-wtaq0*clm%forc_q)

!
! Evaporation flux from foliage
!

     erre = 0.
     if (efe*efeb < 0.0) then
        efeold = efe
        efe  = 0.1*efeold
        erre = efe - efeold
     endif
     dt_veg = (clm%sabv + air + bir*clm%t_veg**4 + cir*tg**4 - efsh - efe) &
          / (- 4.*bir*clm%t_veg**3 +dc1*wtga +dc2*wtgaq*qsatldT)
     clm%t_veg = tlbef + dt_veg

     dels = clm%t_veg-tlbef
     del  = abs(dels)
     err = 0.
     if (del > delmax) then
        dt_veg = delmax*dels/del
        clm%t_veg = tlbef + dt_veg
        err = clm%sabv + air + bir*tlbef**3*(tlbef + 4.*dt_veg) &
             + cir*tg**4 - (efsh + dc1*wtga*dt_veg)          &
             - (efe + dc2*wtgaq*qsatldT*dt_veg)
     endif

!
! Fluxes from leaves to canopy space
! "efe" was limited as its sign changes frequently.  This limit may
! result in an imbalance in "hvap*qflx_evap_veg" and "efe + dc2*wtgaq*qsatldT*dt_veg" 
!

     efpot = clm%forc_rho*wtl*(wtgaq*(qsatl+qsatldT*dt_veg) &
          -wtgq0*qg-wtaq0*clm%forc_q)
     clm%qflx_evap_veg = rpp*efpot

!
! Calculation of evaporative potentials (efpot) and
! interception losses; flux in kg m**-2 s-1.  ecidif 
! holds the excess energy if all intercepted water is evaporated
! during the timestep.  This energy is later added to the
! sensible heat flux.
!

     ecidif = 0.
     if (efpot > 0. .AND. btran > 1.e-10) then
        clm%qflx_tran_veg = efpot*rppdry
     else
        clm%qflx_tran_veg = 0.
     endif
     ecidif = max(0._r4, clm%qflx_evap_veg-clm%qflx_tran_veg-clm%h2ocan/clm%dtime)
     clm%qflx_evap_veg = min(clm%qflx_evap_veg,clm%qflx_tran_veg+clm%h2ocan/clm%dtime)

!
! The energy loss due to above two limits is added to 
! the sensible heat flux.
!

     clm%eflx_sh_veg = efsh + dc1*wtga*dt_veg + err + erre +hvap*ecidif

!
! Re-calculate saturated vapor pressure, specific humidity, and their
! derivatives at the leaf surface
!

     call QSatclm(clm%t_veg, clm%forc_pbot, el, deldT, qsatl, &
               qsatldT    )

!
! Update vegetation/ground surface temperature, canopy air temperature, 
! canopy vapor pressure, aerodynamic temperature, and
! Monin-Obukhov stability parameter for next iteration. 
!

     taf = wtg0*tg + wta0*thm + wtl0*clm%t_veg
     qaf = wtlq0*qsatl+wtgq0*qg+clm%forc_q*wtaq0

!
! Update Monin-Obukhov length and wind speed including the stability effect
!

     dth = thm-taf       
     dqh = clm%forc_q-qaf

     tstar=temp1*dth
     qstar=temp2*dqh

     dthv=dth*(1.+0.61*clm%forc_q)+0.61*th*dqh

     thvstar=tstar*(1.+0.61*clm%forc_q) + 0.61*th*qstar
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
     if (nmozsgn >= 4) then 
        obu = zldis/(-0.01)
     endif

     obuold = obu

!
! Test for convergence
!

     itlef = itlef+1
     if (itlef > itmin) then
        dele = abs(efe-efeb)
        efeb = efe
        det  = max(del,del2)
        if (det < dtmin .AND. dele < dlemin) exit 
     endif

  enddo ITERATION     ! End stability iteration

!
! Energy balance check in canopy
!

  err = clm%sabv + air + bir*tlbef**3*(tlbef + 4.*dt_veg) &
       + cir*tg**4 - clm%eflx_sh_veg - hvap*clm%qflx_evap_veg
  if (abs(err) > 0.1) then
     write(6,*) 'energy balance in canopy X',err
  endif

!
! Fluxes from ground to canopy space 
!

  delt  = wtal*tg-wtl0*clm%t_veg-wta0*thm
  delq  = wtalq*qg-wtlq0*qsatl-wtaq0*clm%forc_q
  clm%taux  = -clm%forc_rho*clm%forc_u/ram(1)
  clm%tauy  = -clm%forc_rho*clm%forc_v/ram(1)
  clm%eflx_sh_grnd = cpair*clm%forc_rho*wtg*delt
  clm%qflx_evap_soi = clm%forc_rho*wtgq*delq

!
! 2 m height air temperature
!

  clm%t_ref2m   = clm%t_ref2m + (taf + temp1*dth * &
       1./vkc *log((2.+z0hv)/z0hv))

!
! Downward longwave radiation below the canopy    
!

  dlrad = (1.-emv)*emg*clm%forc_lwrad &
       + emv*emg * sb * &
       tlbef**3*(tlbef + 4.*dt_veg)

!
! Upward longwave radiation above the canopy    
!

  ulrad = ( (1.-emg)*(1.-emv)*(1.-emv)*clm%forc_lwrad &
       + emv*(1.+(1.-emg)*(1.-emv))*sb * tlbef**3 &
       *(tlbef + 4.*dt_veg) + emg *(1.-emv) *sb * tg**4)

!
! Derivative of soil energy flux with respect to soil temperature (cgrnd) 
!

  cgrnds = cgrnds + cpair*clm%forc_rho*wtg*wtal
  cgrndl = cgrndl + clm%forc_rho*wtgq*wtalq*dqgdT
  cgrnd  = cgrnds  + cgrndl*htvp

!
! Update dew accumulation (kg m-2) 
!

  clm%h2ocan = max(0._r4,clm%h2ocan + (clm%qflx_tran_veg-clm%qflx_evap_veg)*clm%dtime)

end subroutine CanopyFluxes
