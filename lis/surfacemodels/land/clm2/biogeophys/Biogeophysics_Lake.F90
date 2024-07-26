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

subroutine Biogeophysics_Lake (clm) 

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
! Calculates lake temperatures  and surface fluxes.
!
! Method:
! Lake temperatures are determined from a one-dimensional thermal
! stratification model based on eddy diffusion concepts to 
! represent vertical mixing of heat.
!
! d ts    d            d ts     1 ds
! ---- = -- [(km + ke) ----] + -- --
!  dt    dz             dz     cw dz   
!
! where: ts = temperature (kelvin)
!         t = time (s)
!         z = depth (m)
!        km = molecular diffusion coefficient (m**2/s)
!        ke = eddy diffusion coefficient (m**2/s)
!        cw = heat capacity (j/m**3/kelvin)
!         s = heat source term (w/m**2)
!
! There are two types of lakes: 
!   Deep lakes are 50 m. 
!   Shallow lakes are 10 m deep.
!
!   For unfrozen deep lakes:    ke > 0 and    convective mixing
!   For unfrozen shallow lakes: ke = 0 and no convective mixing
!
! Use the Crank-Nicholson method to set up tridiagonal system of equations to
! solve for ts at time n+1, where the temperature equation for layer i is
! r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1
!
! The solution conserves energy as:
!
! cw*([ts(      1)] n+1 - [ts(      1)] n)*dz(      1)/dt + ... +
! cw*([ts(nlevlak)] n+1 - [ts(nlevlak)] n)*dz(nlevlak)/dt = fin
!
! where:
! [ts] n   = old temperature (kelvin)
! [ts] n+1 = new temperature (kelvin)
! fin      = heat flux into lake (w/m**2)
!          = beta*sabg + forc_lwrad - eflx_lwrad_out - eflx_sh_tot - eflx_lh_tot 
!            - hm + phi(1) + ... + phi(nlevlak) 
!
! Author:
! Gordon Bonan
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: Biogeophysics_Lake.F90,v 1.6 2004/11/24 22:56:19 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varpar, only : nlevlak
  use clm2_varcon, only : hvap, hsub, hfus, cpair, cpliq, tkwat, tkice, &
                        sb, vkc, grav, denh2o, tfrz, spval
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm      ! CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer i,j          ! do loop or array index
  integer :: idlak = 1 ! index of lake, 1 = deep lake, 2 = shallow lake
  integer niters       ! maximum number of iterations for surface temperature
  integer iter         ! iteration index
  integer nmozsgn      ! number of times moz changes sign
  real(r8) ax      !
  real(r8) bx      !
  real(r8) beta1   ! coefficient of conective velocity [-]
  real(r8) degdT   ! d(eg)/dT
  real(r8) dqh     ! diff of humidity between ref. height and surface
  real(r8) dth     ! diff of virtual temp. between ref. height and surface
  real(r8) dthv    ! diff of vir. poten. temp. between ref. height and surface
  real(r8) dzsur   !
  real(r8) eg      ! water vapor pressure at temperature T [pa]
  real(r8) emg     ! ground emissivity (0.97 for snow,
  real(r8) hm      ! energy residual [W m-2]
  real(r8) htvp    ! latent heat of vapor of water (or sublimation) [j/kg]
  real(r8) obu     ! monin-obukhov length (m)
  real(r8) obuold  ! monin-obukhov length of previous iteration
  real(r8) qsatg   ! saturated humidity [kg kg-1]
  real(r8) qsatgdT ! d(qsatg)/dT
  real(r8) qstar   ! moisture scaling parameter
  real(r8) ram     ! aerodynamical resistance [s/m]
  real(r8) rah     ! thermal resistance [s/m]
  real(r8) raw     ! moisture resistance [s/m]
  real(r8) stftg3  !
  real(r8) temp1   ! relation for potential temperature profile
  real(r8) temp2   ! relation for specific humidity profile
  real(r8) temp3   ! relation for 2m potential temperature profile
  real(r8) temp4   ! relation for 2m specific humidity profile
  real(r8) tgbef   !
  real(r8) thm     ! intermediate variable (forc_t+0.0098*forc_hgt_t)
  real(r8) thv     ! virtual potential temperature (kelvin)
  real(r8) thvstar ! virtual potential temperature scaling parameter
  real(r8) tksur   ! thermal conductivity of snow/soil (w/m/kelvin)
  real(r8) tstar   ! temperature scaling parameter
  real(r8) um      ! wind speed including the stablity effect [m s-1]
  real(r8) ur      ! wind speed at reference height [m s-1]
  real(r8) ustar   ! friction velocity [m s-1]
  real(r8) wc      ! convective velocity [m s-1]
  real(r8) zeta    ! dimensionless height used in Monin-Obukhov theory
  real(r8) zii     ! convective boundary height [m]
  real(r8) zldis   ! reference height "minus" zero displacement height [m]
  real(r8) displa  ! displacement height [m]
  real(r8) z0mg    ! roughness length over ground, momentum [m]
  real(r8) z0hg    ! roughness length over ground, sensible heat [m]
  real(r8) z0qg    ! roughness length over ground, latent heat [m]
  real(r8) beta(2) ! fraction solar rad absorbed at surface: depends on lake type
  real(r8) za(2)   ! base of surface absorption layer (m): depends on lake type
  real(r8) eta(2)  ! light extinction coefficient (/m): depends on lake type
  real(r8) p0      ! neutral value of turbulent prandtl number
  real(r8) a(nlevlak)    ! "a" vector for tridiagonal matrix
  real(r8) b(nlevlak)    ! "b" vector for tridiagonal matrix
  real(r8) c(nlevlak)    ! "c" vector for tridiagonal matrix
  real(r8) r(nlevlak)    ! "r" vector for tridiagonal solution
  real(r8) rhow(nlevlak) ! density of water (kg/m**3)
  real(r8) phi(nlevlak)  ! solar radiation absorbed by layer (w/m**2)
  real(r8) kme(nlevlak)  ! molecular + eddy diffusion coefficient (m**2/s)
  real(r8) cwat    ! specific heat capacity of water (j/m**3/kelvin)
  real(r8) ws      ! surface friction velocity (m s-1)
  real(r8) ks      ! coefficient
  real(r8) in      ! relative flux of solar radiation into layer
  real(r8) out     ! relative flux of solar radiation out of layer
  real(r8) ri      ! richardson number
  real(r8) fin     ! heat flux into lake - flux out of lake (w/m**2)
  real(r8) ocvts   ! (cwat*(t_lake[n  ])*dz
  real(r8) ncvts   ! (cwat*(t_lake[n+1])*dz
  real(r8) m1      ! intermediate variable for calculating r, a, b, c
  real(r8) m2      ! intermediate variable for calculating r, a, b, c
  real(r8) m3      ! intermediate variable for calculating r, a, b, c
  real(r8) ke      ! eddy diffusion coefficient (m**2/s)
  real(r8) km      ! molecular diffusion coefficient (m**2/s)
  real(r8) zin     ! depth at top of layer (m)
  real(r8) zout    ! depth at bottom of layer (m)
  real(r8) drhodz  ! d [rhow] /dz (kg/m**4)
  real(r8) n2      ! brunt-vaisala frequency (/s**2)
  real(r8) num     ! used in calculating ri
  real(r8) den     ! used in calculating ri
  real(r8) tav     ! used in aver temp for convectively mixed layers
  real(r8) nav     ! used in aver temp for convectively mixed layers
  real(r8) phidum  ! temporary value of phi
  real(r8) u2m     ! 2 m wind speed (m s-1)
  real(r8) cf      ! s m**2/umol -> s/m

!----End Variable List--------------------------------------------------
#if(defined COUPLED)
!!  clm%forc_hgt  =clm%forc_hgt+clm%displa+clm%z0m
  clm%forc_hgt_u=clm%forc_hgt+clm%displa+clm%z0m
  clm%forc_hgt_t=clm%forc_hgt+clm%displa+clm%z0m
  clm%forc_hgt_q=clm%forc_hgt+clm%displa+clm%z0m
#else
  clm%forc_hgt  =10.0+clm%displa+clm%z0m
  clm%forc_hgt_u=10.0+clm%displa+clm%z0m
  clm%forc_hgt_t=2.0+clm%displa+clm%z0m
  clm%forc_hgt_q=2.0+clm%displa+clm%z0m
#endif
!
! Surface Radiation
!

  call SurfaceRadiation (clm)


!
! Determine beginning water balance
!

  clm%begwb = clm%h2osno  

!
! [1] Constants and model parameters
!

!
! Constants for lake temperature model
!

  beta = (/0.4, 0.4/)                            ! (deep lake, shallow lake)
  za   = (/0.6, 0.5/)    
  eta  = (/0.1, 0.5/)  
  p0   = 1.  

!
! Roughness lengths
!

  if (clm%t_grnd >= tfrz)then                    ! for unfrozen lake
     z0mg = 0.01
  else                                           ! for frozen lake
     z0mg = 0.04
  endif
  z0hg = z0mg
  z0qg = z0mg

!
! Latent heat 
!

  if (clm%forc_t > tfrz) then
     htvp = hvap
  else
     htvp = hsub
  endif

#if (defined PERGRO)
  htvp = hvap
#endif

!
! Emissivity
!

  emg = 0.97

!
! [2] Surface temperature and fluxes
!

  dzsur = clm%dz(1) + clm%snowdp

!
! Saturated vapor pressure, specific humidity and their derivatives
! at lake surface
!

  call QSatclm(clm%t_grnd, clm%forc_pbot, eg, degdT, qsatg, &
            qsatgdT     )

!
! Potential, virtual potential temperature, and wind speed at the
! reference height
!

  beta1=1.       ! -  (in computing W_*)
  zii = 1000.    ! m  (pbl height)
  thm = clm%forc_t + 0.0098*clm%forc_hgt_t         ! intermediate variable 
  thv = clm%forc_th*(1.+0.61*clm%forc_q)           ! virtual potential T
  ur = max(1.0_r4,sqrt(clm%forc_u*clm%forc_u+clm%forc_v*clm%forc_v))

!
! Initialize stability variables
!

  nmozsgn = 0
  obuold = 0.

  dth   = thm-clm%t_grnd
  dqh   = clm%forc_q-qsatg
  dthv  = dth*(1.+0.61*clm%forc_q)+0.61*clm%forc_th*dqh
  zldis = clm%forc_hgt_u-0.

!
! Initialize Monin-Obukhov length and wind speed
!

  call MoninObukIni(ur, thv, dthv, zldis, z0mg, &
                    um, obu  )

!
! Begin stability iteration
!

  niters = 3
  do iter = 1, niters
     tgbef = clm%t_grnd
     if (clm%t_grnd > tfrz) then
        tksur = tkwat
     else
        tksur = tkice
     endif

!
! Determine friction velocity, and potential temperature and humidity
! profiles of the surface boundary layer
!

     displa = 0.0_r4
     call FrictionVelocity (displa, z0mg,  z0hg,  z0qg,  obu, &
                            iter, ur, um, ustar, temp1, temp2, &
                            temp3,temp4,clm)
     obuold = obu

!
! Determine aerodynamic resistances
!
      clm%ch = ustar*ustar/um ! Add-in for ALMA output
      clm%chs2 = temp3*ustar ! For calc of 2M T
      clm%cqs2 = temp4*ustar ! For calc of 2M Q
 
     ram    = 1./(ustar*ustar/um)
     rah    = 1./(temp1*ustar)
     raw    = 1./(temp2*ustar)

!
! Get derivative of fluxes with respect to ground temperature
!

     stftg3 = emg*sb*tgbef*tgbef*tgbef

     ax  = clm%sabg + emg*clm%forc_lwrad + 3.*stftg3*tgbef &
          + clm%forc_rho*cpair/rah*thm &
          - htvp*clm%forc_rho/raw*(qsatg-qsatgdT*tgbef - clm%forc_q) &
          + tksur*clm%t_lake(1)/dzsur

     bx  = 4.*stftg3 + clm%forc_rho*cpair/rah &
          + htvp*clm%forc_rho/raw*qsatgdT + tksur/dzsur

     clm%t_grnd = ax/bx

!
! Surface fluxes of momentum, sensible and latent heat
! using ground temperatures from previous time step
!

     clm%eflx_sh_grnd = clm%forc_rho*cpair*(clm%t_grnd-thm)/rah
     clm%qflx_evap_soi = clm%forc_rho*(qsatg+qsatgdT*(clm%t_grnd-tgbef)-clm%forc_q)/raw

!
! Re-calculate saturated vapor pressure, specific humidity and their
! derivatives at lake surface
!

     call QSatclm(clm%t_grnd, clm%forc_pbot, eg, degdT, qsatg, &
               qsatgdT     )

     dth=thm-clm%t_grnd
     dqh=clm%forc_q-qsatg

     tstar = temp1*dth
     qstar = temp2*dqh

     dthv=dth*(1.+0.61*clm%forc_q)+0.61*clm%forc_th*dqh
     thvstar=tstar*(1.+0.61*clm%forc_q) + 0.61*clm%forc_th*qstar
     zeta=zldis*vkc * grav*thvstar/(ustar**2*thv)

     if (zeta >= 0.) then     !stable
        zeta = min(2._r4,max(zeta,0.01_r4))
        um = max(ur,0.1_r4)
     else                     !unstable
        zeta = max(-100._r4,min(zeta,-0.01_r4))
        wc = beta1*(-grav*ustar*thvstar*zii/thv)**0.333
        um = sqrt(ur*ur+wc*wc)
     endif
     obu = zldis/zeta

     if (obuold*obu < 0.) nmozsgn = nmozsgn+1
     if (nmozsgn >= 4) EXIT

  enddo

!
! If there is snow on the ground and t_grnd > tfrz: reset t_grnd = tfrz.
! Re-evaluate ground fluxes. Energy inbalance used to melt snow.  
! h2osno > 0.5 prevents spurious fluxes.
!

  if (clm%h2osno > 0.5 .AND. clm%t_grnd > tfrz) then
     clm%t_grnd = tfrz
     clm%eflx_sh_grnd = clm%forc_rho*cpair*(clm%t_grnd-thm)/rah
     clm%qflx_evap_soi = clm%forc_rho*(qsatg+qsatgdT*(clm%t_grnd-tgbef) & 
                       - clm%forc_q)/raw  !note that qsatg and qsatgdT should be f(tgbef)
  endif

!
! Net longwave from ground to atmosphere
!

  clm%eflx_lwrad_out = (1.-emg)*clm%forc_lwrad + stftg3*(-3.*tgbef+4.*clm%t_grnd)

!
! Radiative temperature
!

  clm%t_rad = (clm%eflx_lwrad_out/sb)**0.25

!
! Ground heat flux
!

  clm%eflx_soil_grnd = clm%sabg + clm%forc_lwrad - clm%eflx_lwrad_out - &
                       clm%eflx_sh_grnd - htvp*clm%qflx_evap_soi

  clm%taux   = -clm%forc_rho*clm%forc_u/ram
  clm%tauy   = -clm%forc_rho*clm%forc_v/ram

  clm%eflx_sh_tot   = clm%eflx_sh_grnd
  clm%qflx_evap_tot = clm%qflx_evap_soi
  clm%eflx_lh_tot   = htvp*clm%qflx_evap_soi
  clm%eflx_lh_grnd  = htvp*clm%qflx_evap_soi

!
! 2 m height air temperature
!

  clm%t_ref2m   = (clm%t_grnd + temp1*dth * 1./ &
                  vkc *log((2.+z0hg)/z0hg))

!
! Energy residual used for melting snow
!

  if (clm%h2osno > 0. .AND. clm%t_grnd >= tfrz) then
     hm = min( clm%h2osno*hfus/clm%dtime, max(clm%eflx_soil_grnd,0._r4) )
  else
     hm = 0.
  endif
  clm%qmelt = hm/hfus             ! snow melt (mm s-1)

!
! [3] Lake layer temperature
!

!
! Lake density
!

  do j = 1, nlevlak
     rhow(j) = 1000.*( 1.0 - 1.9549e-05*(abs(clm%t_lake(j)-277.))**1.68 )
  enddo

!
! Eddy diffusion +  molecular diffusion coefficient:
! eddy diffusion coefficient used for unfrozen deep lakes only
!

  cwat = cpliq*denh2o
  km = tkwat/cwat

  fin = beta(idlak) * clm%sabg + clm%forc_lwrad - (clm%eflx_lwrad_out + &
        clm%eflx_sh_tot + clm%eflx_lh_tot + hm)
  u2m = max(1.0_r4,ustar/vkc*log(2./z0mg))

  ws = 1.2e-03 * u2m
  ks = 6.6*sqrt(abs(sin(clm%lat)))*(u2m**(-1.84))

  do j = 1, nlevlak-1
     drhodz = (rhow(j+1)-rhow(j)) / (clm%z(j+1)-clm%z(j))
     n2 = -grav / rhow(j) * drhodz
     num = 40. * n2 * (vkc*clm%z(j))**2
     den = max( (ws**2) * exp(-2.*ks*clm%z(j)), 1.e-10_r4 )
     ri = ( -1. + sqrt( max(1.+num/den, 0._r4) ) ) / 20.
     if (idlak == 1 .AND. clm%t_grnd > tfrz) then
        ke = vkc*ws*clm%z(j)/p0 * exp(-ks*clm%z(j)) / (1.+37.*ri*ri)
     else
        ke = 0.
     endif
     kme(j) = km + ke 
  enddo

  kme(nlevlak) = kme(nlevlak-1)

!
! Heat source term: unfrozen lakes only
!

  do j = 1, nlevlak
     zin  = clm%z(j) - 0.5*clm%dz(j)
     zout = clm%z(j) + 0.5*clm%dz(j)
     in  = exp( -eta(idlak)*max(  zin-za(idlak),0._r4 ) )
     out = exp( -eta(idlak)*max( zout-za(idlak),0._r4 ) )

!
! Assume solar absorption is only in the considered depth
!

     if (j == nlevlak) out = 0.  
     if (clm%t_grnd > tfrz) then
        phidum = (in-out) * clm%sabg * (1.-beta(idlak))
     else if (j == 1) then
        phidum= clm%sabg * (1.-beta(idlak))
     else
        phidum = 0.
     endif
     phi(j) = phidum
  enddo

!
! Sum cwat*t_lake*dz for energy check
!

  ocvts = 0.
  do j = 1, nlevlak
     ocvts = ocvts + cwat*clm%t_lake(j)*clm%dz(j) 
  enddo

!
! Set up vector r and vectors a, b, c that define tridiagonal matrix
!

  j = 1
  m2 = clm%dz(j)/kme(j) + clm%dz(j+1)/kme(j+1)
  m3 = clm%dtime/clm%dz(j)
  r(j) = clm%t_lake(j) + (fin+phi(j))*m3/cwat - (clm%t_lake(j)-clm%t_lake(j+1))*m3/m2
  a(j) = 0.
  b(j) = 1. + m3/m2
  c(j) = -m3/m2

  j = nlevlak
  m1 = clm%dz(j-1)/kme(j-1) + clm%dz(j)/kme(j)
  m3 = clm%dtime/clm%dz(j)
  r(j) = clm%t_lake(j) + phi(j)*m3/cwat + (clm%t_lake(j-1)-clm%t_lake(j))*m3/m1
  a(j) = -m3/m1
  b(j) = 1. + m3/m1
  c(j) = 0.

  do j = 2, nlevlak-1
     m1 = clm%dz(j-1)/kme(j-1) + clm%dz(j  )/kme(j  )
     m2 = clm%dz(j  )/kme(j  ) + clm%dz(j+1)/kme(j+1)
     m3 = clm%dtime/clm%dz(j)
     r(j) = clm%t_lake(j) + phi(j)*m3/cwat + &
            (clm%t_lake(j-1) - clm%t_lake(j))*m3/m1 - &
            (clm%t_lake(j)-clm%t_lake(j+1))*m3/m2

     a(j) = -m3/m1
     b(j) = 1. + m3/m1 + m3/m2
     c(j) = -m3/m2
  enddo

!
! Solve for t_lake: a, b, c, r, u 
!

  call Tridiagonal (nlevlak, a, b, c, r, clm%t_lake(1:nlevlak)) 

!
! Convective mixing: make sure cwat*dz*ts is conserved.
!

  if (idlak == 1 .AND. clm%t_grnd > tfrz) then
     do j = 1, nlevlak-1
        if (rhow(j) > rhow(j+1)) then
           tav = 0.
           nav = 0.
           do i = 1, j+1
              tav = tav + clm%t_lake(i)*clm%dz(i)
              nav = nav + clm%dz(i)
           enddo
           tav = tav/nav
           do i = 1, j+1
              clm%t_lake(i) = tav
              rhow(i) = 1000.*( 1.0 - 1.9549e-05*(abs(clm%t_lake(i)-277.))**1.68 )
           enddo
        endif
     enddo
  endif
 
!
! Sum cwat*t_lake*dz and total energy into lake for energy check
!

  ncvts = 0.
  do j = 1, nlevlak
     ncvts = ncvts + cwat*clm%t_lake(j)*clm%dz(j) 
     fin = fin + phi(j)
  enddo

!  clm%errsoi = (ncvts-ocvts) / clm%dtime - fin

!
! [4] Set other clm values for lake points
!

!
! The following are needed for global average on history tape.
! Note: time invariant variables set in initialization phase:
! z, dz, snl, h2osoi_liq, and h2osoi_ice
!

  clm%t_veg = clm%forc_t  ! to be consistent with treatment of t_veg for bare soil points

  clm%eflx_sh_veg     = 0.
  clm%eflx_lh_vegt    = 0.
  clm%eflx_lh_vege    = 0.
  clm%eflx_lwrad_net  = clm%eflx_lwrad_out -  clm%forc_lwrad
  
! Components that are not displayed over lake on history tape and 
! therefore need to be set to spval here

  clm%rssun    = spval
  clm%rssha    = spval
  clm%t_snow   = spval

end subroutine Biogeophysics_Lake
