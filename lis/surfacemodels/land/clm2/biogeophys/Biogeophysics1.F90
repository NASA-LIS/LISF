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
subroutine Biogeophysics1 (clm,cgrnd,cgrndl,cgrnds,tg,emg,htvp, dlrad,ulrad,tssbef)

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
! This is the main subroutine to execute the calculation of leaf temperature
! and surface fluxes. Biogeophysics2.F90 then determines soil/snow and ground
! temperatures and updates the surface fluxes for the new ground
! temperature.
!
! Calling sequence is:
!  Biogeophysics1:                   surface biogeophysics driver
!    -> QSatclm:                        saturated vapor pressure, specific humidity,
!                                     and derivatives at ground surface
!    -> SurfaceRadiation:            surface solar radiation
!    -> BareGroundFluxes:            surface fluxes for bare soil or
!                                     snow-covered vegetation patches
!          -> MoninObukIni:          first-guess Monin-Obukhov length and
!                                     wind speed
!          -> FrictionVelocity:      friction velocity and potential 
!                                     temperature and humidity profiles
!    -> CanopyFluxes:                leaf temperature and surface fluxes
!                                     for vegetated patches 
!          -> QSatclm                   saturated vapor pressure, specific humidity,
!                                     and derivatives at leaf surface
!          -> MoninObukIni           first-guess Monin-Obukhov length and 
!                                     wind speed
!          -> FrictionVelocity       friction velocity and potential
!                                     temperature and humidity profiles
!          -> Stomata                stomatal resistance and photosynthesis
!                                     for sunlit leaves
!          -> Stomata                stomatal resistance and photosynthesis
!                                     for shaded leaves
!          -> SensibleHCond          sensible heat conductance for air, leaf,
!                                     and ground
!          -> LatentHCond            latent heat conductance for ground and
!                                     leaf
!          -> QSatclm                   recalculation of saturated vapor pressure,
!                                     specific humidity, and derivatives at
!                                     leaf surface using updated leaf temperature
!  Leaf temperature
!   Foliage energy conservation is given by the foliage energy budget 
!   equation:
!                  Rnet - Hf - LEf = 0 
!   The equation is solved by Newton-Raphson iteration, in which this 
!   iteration includes the calculation of the photosynthesis and 
!   stomatal resistance, and the integration of turbulent flux profiles. 
!   The sensible and latent heat transfer between foliage and atmosphere 
!   and ground is linked by the equations:  
!                  Ha = Hf + Hg and Ea = Ef + Eg
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: Biogeophysics1.F90,v 1.8 2004/11/24 22:56:17 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, only : denh2o, denice, roverg, hvap, hsub, istice, istwet 
  use clm2_varpar, only : nlevsoi 

  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout)  :: clm !CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer i      ! loop indices
  real(r8) qred    ! soil surface relative humidity
  real(r8) avmuir  ! ir inverse optical depth per unit leaf area
  real(r8) eg      ! water vapor pressure at temperature T [pa]
  real(r8) qsatg   ! saturated humidity [kg kg-1]
  real(r8) degdT   ! d(eg)/dT
  real(r8) qsatgdT ! d(qsatg)/dT
  real(r8) fac     ! soil wetness of surface layer
  real(r8) psit    ! negative potential of soil
  real(r8) hr      ! relative humidity
  real(r8) wx      ! partial volume of ice and water of surface layer
  real(r8) zlnd
  real(r8) zsno
  real(r8) smpmin
  real(r8) qg
  real(r8) dqgdT
  real(r8) emv
  real(r8) z0mg
  real(r8) z0hg
  real(r8) z0qg
  real(r8) z0mv
  real(r8) z0hv
  real(r8) z0qv
  real(r8) beta
  real(r8) zii
  real(r8) thm
  real(r8) thv
  real(r8) ur

  real(r8) :: cgrnd
  real(r8) :: cgrndl
  real(r8) :: cgrnds
  real(r8) :: tg
  real(r8) :: emg
  real(r8) :: htvp
  real(r8) :: dlrad
  real(r8) :: ulrad
  real(r8) :: tssbef(-5:10)
!----End Variable List--------------------------------------------------

!
! Initial set 
!
  zlnd = 0.01
  zsno = 0.0024
  smpmin = -1.e8
  clm%eflx_sh_tot    = 0.
  clm%qflx_evap_tot  = 0.
  clm%eflx_lh_tot    = 0.
  clm%eflx_sh_veg    = 0.  
  clm%qflx_evap_veg  = 0.  
  clm%qflx_tran_veg  = 0.  
  cgrnd          = 0._r8
  cgrnds         = 0._r8
  cgrndl         = 0._r8
  clm%t_ref2m        = 0.
!
! Ground and soil temperatures from previous time step
!
  tg = clm%t_soisno(clm%snl+1)
  do i = clm%snl+1, nlevsoi
     tssbef(i) = clm%t_soisno(i)
  enddo

!
! Saturated vapor pressure, specific humidity and their derivatives
! at ground surface
!
  qred = 1.
  if (clm%itypwat/=istwet .AND. clm%itypwat/=istice) then
     wx   = (clm%h2osoi_liq(1)/denh2o+clm%h2osoi_ice(1)/denice)/clm%dz(1)
     fac  = min(1._r4, wx/clm%watsat(1))
     fac  = max( fac, 0.01_r4 )
     psit = -clm%sucsat(1) * fac ** (- clm%bsw(1))
     psit = max(smpmin, psit)
     hr   = exp(psit/roverg/tg)
     qred = (1.-clm%frac_sno)*hr + clm%frac_sno
  endif

  call QSatclm(tg, clm%forc_pbot, eg, degdT, qsatg, &
            qsatgdT)
  qg = qred*qsatg  
  dqgdT = qred*qsatgdT

  if (qsatg > clm%forc_q .AND. clm%forc_q > qred*qsatg) then
     qg = clm%forc_q
     dqgdT = 0.
  endif

!
! Emissivity
!

  if (clm%h2osno>0. .OR.clm%itypwat==istice) then
     emg = 0.97
  else
     emg = 0.96
  endif
  avmuir=1.
  emv=1.-exp(-(clm%elai+clm%esai)/avmuir)

!
! Latent heat. We arbitrarily assume that the sublimation occurs 
! only as h2osoi_liq = 0
!

  htvp = hvap
  if (clm%h2osoi_liq(clm%snl+1) <= 0. .AND. clm%h2osoi_ice(clm%snl+1) > 0.) htvp = hsub

!
! Switch between vaporization and sublimation causes rapid solution
! separation in perturbation growth test
!

#if (defined PERGRO)
  htvp = hvap
#endif

!
! Roughness lengths
!

  if (clm%frac_sno > 0.) then
     z0mg = zsno
     z0hg = z0mg            ! initial set only
     z0qg = z0mg            ! initial set only
  else
     z0mg = zlnd
     z0hg = z0mg            ! initial set only
     z0qg = z0mg            ! initial set only
  endif

  clm%z0m = clm%z0mr*clm%htop
  clm%displa = clm%displar*clm%htop
!!print*, 'displa, displar, htop:', clm%displa, clm%displar, clm%htop 
  z0mv = clm%z0m
  z0hv = z0mv
  z0qv = z0mv


!=== LDAS modifications****
!=== New code added for LDAS framework
!=== Added specification of forcing heights here instead of in the ldasdrv.f
!=== since displacement heights set here as opposed to being read in from a
!=== text file

!=== This change was necessary because the forcing heights would be less than the
!=== canopy heights and cause model crash 
#if ( defined COUPLED )
!!clm%forc_hgt  =clm%forc_hgt+clm%displa+clm%z0m
  clm%forc_hgt_u=clm%forc_hgt+clm%displa+clm%z0m
  clm%forc_hgt_t=clm%forc_hgt+clm%displa+clm%z0m
  clm%forc_hgt_q=clm%forc_hgt+clm%displa+clm%z0m
!!print*, 'forc_hgt', clm%forc_hgt, clm%forc_hgt_u, clm%forc_hgt_t, clm%forc_hgt_q
#else
  clm%forc_hgt  =10.0+clm%displa+clm%z0m
  clm%forc_hgt_u=10.0+clm%displa+clm%z0m
  clm%forc_hgt_t=2.0+clm%displa+clm%z0m
  clm%forc_hgt_q=2.0+clm%displa+clm%z0m
#endif
!    print*, 'in bio..',clm%t_veg

!
! Potential, virtual potential temperature, and wind speed at the 
! reference height
!

  beta=1.
  zii = 1000.
  thm = clm%forc_t + 0.0098*clm%forc_hgt_t              
  thv = clm%forc_th*(1.+0.61*clm%forc_q)
  ur = max(1.0_r8,sqrt(clm%forc_u*clm%forc_u+clm%forc_v*clm%forc_v))
!  print*, ur

!
! Surface Radiation
!

  call SurfaceRadiation (clm)
!
! Surface Temperature and Fluxes
!

!
! BARE SOIL OR SNOW-COVERED VEGETATION
! Ground fluxes
! NOTE: in the current scheme clm%frac_veg_nosno is EITHER 1 or 0
!

!  print*, clm%itypveg,clm%t_veg,clm%t_grnd,clm%displa,clm%z0m,
!  write(*,33) clm%itypveg,clm%t_veg,clm%t_grnd,clm%displa,clm%z0m,thm,&
!& clm%forc_th,thm,clm%tg,clm%qg,clm%dqgdT,emv,clm%emg,clm%dlrad,clm%ulrad,&
!& clm%cgrnds,clm%cgrndl,clm%cgrnd
! 33 format(i3,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,3(f8.4,1x),3(f8.4,1x),1x,7(f8.4,1x))

  if (clm%frac_veg_nosno == 0) then

     call BareGroundFluxes (tg,     thm,   qg,    thv,   z0mg,   &
                            z0hg,   z0qg,  dqgdT, htvp,  beta,   &
                            zii,    ur,    dlrad, ulrad, cgrnds, &
                            cgrndl, cgrnd, clm    )
     clm%psnsun = 0.
     clm%psnsha = 0. !put these lines here to avoid psn = NaN

!         print*, 'in bio..',clm%t_veg

!
! VEGETATION
! Calculate canopy temperature, latent and sensible fluxes from the canopy,
! and leaf water change by evapotranspiration
!

  else

     call CanopyFluxes (z0mv,   z0hv,  z0qv,  thm,   clm%forc_th, &
                        thv,    tg,    qg,    dqgdT, htvp,        &
                        emv,    emg,   dlrad, ulrad, cgrnds,      &
                        cgrndl, cgrnd, clm    )

!         print*, clm%t_veg, clm%t_grnd, clm%displa, clm%z0m, clm%frac_veg_nosno

  endif
  
!  print*, clm%itypveg, clm%t_veg, clm%t_grnd, clm%displa, clm%z0m, clm%frac_veg_nosno
  

end subroutine Biogeophysics1
