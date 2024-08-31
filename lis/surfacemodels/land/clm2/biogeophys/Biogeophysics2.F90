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

subroutine Biogeophysics2 (clm,cgrnd,cgrndl,cgrnds,tg,emg,htvp, dlrad,ulrad,tssbef)

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
! This is the main subroutine to execute the calculation of soil/snow and
! ground temperatures and update surface fluxes based on the new ground
! temperature 
! 
! Method:
! Calling sequence is:
! Biogeophysics2:                    surface biogeophysics driver
!    -> SoilTemperature:             soil/snow and ground temperatures      
!          -> SoilTermProp           thermal conductivities and heat 
!                                     capacities        
!          -> Tridiagonal            tridiagonal matrix solution            
!          -> PhaseChange            phase change of liquid/ice contents        
!
! (1) Snow and soil temperatures
!     o The volumetric heat capacity is calculated as a linear combination 
!       in terms of the volumetric fraction of the constituent phases. 
!     o The thermal conductivity of soil is computed from 
!       the algorithm of Johansen (as reported by Farouki 1981), and the 
!       conductivity of snow is from the formulation used in
!       SNTHERM (Jordan 1991).
!     o Boundary conditions:  
!       F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
!     o Soil / snow temperature is predicted from heat conduction 
!       in 10 soil layers and up to 5 snow layers. 
!       The thermal conductivities at the interfaces between two 
!       neighboring layers (j, j+1) are derived from an assumption that 
!       the flux across the interface is equal to that from the node j 
!       to the interface and the flux from the interface to the node j+1. 
!       The equation is solved using the Crank-Nicholson method and 
!       results in a tridiagonal system equation.
!
! (2) Phase change (see PhaseChange.F90)
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: Biogeophysics2.F90,v 1.8 2004/11/24 22:56:18 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, only : hvap, cpair, grav, vkc, tfrz, sb
  use clm2_varpar, only : nlevsoi
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout)  :: clm !CLM 1-D Module

!----Local Variables----------------------------------------------------

!  integer j                          ! do loop index
  real(r8) fact(clm%snl+1 : nlevsoi) ! used in computing tridiagonal matrix
  real(r8) egsmax ! max. evaporation which soil can provide at one time step
  real(r8) egidif ! the excess of evaporation over "egsmax"
  real(r8) xmf    ! total latent heat of phase change of ground water
  real(r8) tinc   ! temperature difference of two time step

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
! Determine soil temperatures including surface soil temperature
!

  call SoilTemperature(clm      , tssbef, htvp, emg, cgrnd, &
                       dlrad, tg    , xmf     , fact )

!
! Correct fluxes to present soil temperature
!

  tinc = clm%t_soisno(clm%snl+1) - tssbef(clm%snl+1)
  clm%eflx_sh_grnd =  clm%eflx_sh_grnd + tinc*cgrnds 
  clm%qflx_evap_soi =  clm%qflx_evap_soi + tinc*cgrndl

!
! egidif holds the excess energy if all water is evaporated from
! the top soil layer during the timestep.  This energy is added to
! the sensible heat flux.
!

  egsmax = (clm%h2osoi_ice(clm%snl+1)+clm%h2osoi_liq(clm%snl+1)) / clm%dtime
  egidif = max( 0._r4, clm%qflx_evap_soi - egsmax )
  clm%qflx_evap_soi = min ( clm%qflx_evap_soi, egsmax )
  clm%eflx_sh_grnd = clm%eflx_sh_grnd + htvp*egidif

!
! Ground heat flux
!

  clm%eflx_soil_grnd = clm%sabg + dlrad + (1-clm%frac_veg_nosno)*emg*clm%forc_lwrad &
       - emg*sb*tssbef(clm%snl+1)**3*(tssbef(clm%snl+1) + 4.*tinc) &
       - (clm%eflx_sh_grnd+clm%qflx_evap_soi*htvp)

!
! Total fluxes (vegetation + ground)
!
  clm%eflx_sh_tot = clm%eflx_sh_veg + clm%eflx_sh_grnd
  clm%qflx_evap_tot = clm%qflx_evap_veg + clm%qflx_evap_soi
  clm%eflx_lh_tot= hvap*clm%qflx_evap_veg + htvp*clm%qflx_evap_soi   ! (account for sublimation)


!
! Assign ground evaporation to sublimation from soil ice or to dew
! on snow or ground 
!

  clm%qflx_evap_grnd = 0.
  clm%qflx_sub_snow = 0.
  clm%qflx_dew_snow = 0.
  clm%qflx_dew_grnd = 0.

  if (clm%qflx_evap_soi >= 0.) then
     ! Do not allow for sublimation in melting (melting ==> evap. ==> sublimation)
     clm%qflx_evap_grnd = min(clm%h2osoi_liq(clm%snl+1)/clm%dtime, clm%qflx_evap_soi)
     clm%qflx_sub_snow = clm%qflx_evap_soi - clm%qflx_evap_grnd
  else
     if (tg < tfrz) then
        clm%qflx_dew_snow = abs(clm%qflx_evap_soi)
     else
        clm%qflx_dew_grnd = abs(clm%qflx_evap_soi)
     endif
  endif

!
! Outgoing long-wave radiation from vegetation + ground
!

  clm%eflx_lwrad_out = ulrad &
       + (1-clm%frac_veg_nosno)*(1.-emg)*clm%forc_lwrad &
       + (1-clm%frac_veg_nosno)*emg*sb * tssbef(clm%snl+1)**4 &
       ! For conservation we put the increase of ground longwave to outgoing
       + 4.*emg*sb*tssbef(clm%snl+1)**3*tinc


!
! Radiative temperature
!

  clm%t_rad = (clm%eflx_lwrad_out/sb)**0.25

!
! Soil Energy balance check
!

!  clm%errsoi = 0. 
!  do j = clm%snl+1, nlevsoi
!     clm%errsoi = clm%errsoi - (clm%t_soisno(j)-tssbef(j))/fact(j) 
!  enddo
!  clm%errsoi = clm%errsoi + clm%eflx_soil_grnd - xmf

!
! Variables needed by history tape
!

! clm%dt_grnd        = tinc
! clm%eflx_lh_vege   = (clm%qflx_evap_veg - clm%qflx_tran_veg) * hvap
! clm%eflx_lh_vegt   = clm%qflx_tran_veg * hvap       
! clm%eflx_lh_grnd   = clm%qflx_evap_soi * htvp
 clm%eflx_lwrad_net = clm%eflx_lwrad_out -  clm%forc_lwrad  

end subroutine Biogeophysics2
