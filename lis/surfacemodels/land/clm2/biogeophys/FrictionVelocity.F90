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

subroutine FrictionVelocity (displa, z0m,   z0h,   z0q,   obu, &
                             iter, ur, um, ustar, temp1, temp2, &
                             temp3, temp4, clm)

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
! Calculation of the friction velocity, relation for potential 
! temperature and humidity profiles of surface boundary layer. 
!
! Method:
! The scheme is based on the work of Zeng et al. (1998): 
! Intercomparison of bulk aerodynamic algorithms for the computation 
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, 
! Vol. 11, 2628-2644.
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: FrictionVelocity.F90,v 1.6 2004/11/24 22:56:25 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, only : vkc
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm !CLM 1-D Module

  real(r8), intent(in) :: displa ! displacement height [m]
  real(r8), intent(in) :: z0m    ! roughness length, momentum [m]
  real(r8), intent(in) :: z0h    ! roughness length, sensible heat [m]
  real(r8), intent(in) :: z0q    ! roughness length, latent heat [m]
  real(r8), intent(in) :: obu    ! monin-obukhov length (m)
  real(r8), intent(in) :: um     ! wind speed including the stablity effect [m s-1]
  real(r8), intent(in) :: ur
  integer , intent(in) :: iter

  real(r8), intent(out) :: ustar ! friction velocity [m s-1]
  real(r8), intent(out) :: temp1 ! relation for potential temperature profile
  real(r8), intent(out) :: temp2 ! relation for specific humidity profile
  real(r8), intent(out) :: temp3 ! relation for 2m potential temperature profile
  real(r8), intent(out) :: temp4 ! relation for 2m specific humidity profile
 
!----Local Variables----------------------------------------------------

  real(r8) zldis   ! reference height "minus" zero displacement heght [m]
  real(r8) StabilityFunc ! stability function for unstable case
  real(r8) zetam   ! transition point of flux-gradient relation (wind profile)
  real(r8) zetat   ! transition point of flux-gradient relation (temp. profile)
  real(r8) zeta    ! dimensionless height used in Monin-Obukhov theory
#if (defined BGC)
! Variables used to diagnose the 10 meter wind
  real(r8) :: tmp1,tmp2,tmp3,tmp4
  real(r8) :: fmnew
  real(r8) :: fm10
  real(r8) :: zeta10
  real(r8), save :: fm
#endif

!----End Variable List--------------------------------------------------

!
! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.
! Wind profile
!

  zldis=clm%forc_hgt_u-displa
  zeta=zldis/obu
  zetam=1.574

  if (zeta < -zetam) then           ! zeta < -1
     ustar=vkc*um/(log(-zetam*obu/z0m)- &
          StabilityFunc(1,-zetam) +StabilityFunc(1,z0m/obu) &
          +1.14*((-zeta)**0.333-(zetam)**0.333))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     ustar=vkc*um/(log(zldis/z0m)- &
          StabilityFunc(1,zeta)+StabilityFunc(1,z0m/obu))
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     ustar=vkc*um/(log(zldis/z0m) + &
          5.*zeta -5.*z0m/obu)
  else                             !  1 < zeta, phi=5+zeta
     ustar=vkc*um/(log(obu/z0m)+5.-5.*z0m/obu &
          +(5.*log(zeta)+zeta-1.))
  endif

#if (defined PERGRO)
  if (zeta < -zetam) then           ! zeta < -1
     ustar=vkc*um/log(-zetam*obu/z0m)
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     ustar=vkc*um/log(zldis/z0m)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     ustar=vkc*um/log(zldis/z0m)
  else                             !  1 < zeta, phi=5+zeta
     ustar=vkc*um/log(obu/z0m)
  endif
#endif


#if (defined BGC)
!
! diagnose 10-m wind for dust model (dstmbl.F)
! Notes from C. Zender's dst.F:
! According to Bon96 p. 62, the displacement height d (here displa) is
! 0.0 <= d <= 0.34 m in dust source regions (i.e., regions w/o trees).
! Therefore d <= 0.034*z1 and may safely be neglected.
! Code from LSM routine SurfaceTemperature was used to obtain u10
!
  if (min(zeta,1.) < 0.) then
     tmp1 = (1. - 16.*min(zeta,1.))**0.25
     tmp2 = log((1.+tmp1*tmp1)/2.)
     tmp3 = log((1.+tmp1)/2.)
     fmnew = 2.*tmp3 + tmp2 - 2.*atan(tmp1) + 1.5707963
  else
     fmnew = -5.*min(zeta,1.)
  endif
  if (iter == 1) then
     fm = fmnew
  else
     fm = 0.5 * (fm+fmnew)
  endif

  zeta10 = min(10./obu, 1.)
  if (zeta == 0.) zeta10 = 0.

  if (zeta10 < 0.) then
     tmp1 = (1.0 - 16.0 * zeta10)**0.25
     tmp2 = log((1.0 + tmp1*tmp1)/2.0)
     tmp3 = log((1.0 + tmp1)/2.0)
     fm10 = 2.0*tmp3 + tmp2 - 2.0*atan(tmp1) + 1.5707963
  else  ! not stable
     fm10 = -5.0 * zeta10
  endif ! not stable
  tmp4 = log(clm%forc_hgt / 10.)

!  clm%u10 = ur - ustar/vkc * (tmp4 - fm + fm10)
!  clm%fv  = ustar
#endif

!
! Temperature profile
!

  zldis=clm%forc_hgt_t-displa
  zeta=zldis/obu
  zetat=0.465
  if (zeta < -zetat) then           ! zeta < -1
     temp1=vkc/(log(-zetat*obu/z0h)-StabilityFunc(2,-zetat) &
          + StabilityFunc(2,z0h/obu) &
          + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp1=vkc/(log(zldis/z0h) - StabilityFunc(2,zeta) + &
          StabilityFunc(2,z0h/obu))
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp1=vkc/(log(zldis/z0h) + 5.*zeta - 5.*z0h/obu)
  else                             !  1 < zeta, phi=5+zeta
     temp1=vkc/(log(obu/z0h) + 5. - 5.*z0h/obu &
          + (5.*log(zeta)+zeta-1.))
  endif

  ! 2M Calculation
  zldis=2.0 + z0h
  zeta=zldis/obu
  zetat=0.465
  if (zeta < -zetat) then           ! zeta < -1
     temp3=vkc/(log(-zetat*obu/z0h)-StabilityFunc(2,-zetat) &
          + StabilityFunc(2,z0h/obu) &
          + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp3=vkc/(log(zldis/z0h) - StabilityFunc(2,zeta) + &
          StabilityFunc(2,z0h/obu))
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp3=vkc/(log(zldis/z0h) + 5.*zeta - 5.*z0h/obu)
  else                             !  1 < zeta, phi=5+zeta
     temp3=vkc/(log(obu/z0h) + 5. - 5.*z0h/obu &
          + (5.*log(zeta)+zeta-1.))
  endif

#if (defined PERGRO)
  if (zeta < -zetat) then           ! zeta < -1
     temp1=vkc/log(-zetat*obu/z0h)
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp1=vkc/log(zldis/z0h)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp1=vkc/log(zldis/z0h)
  else                             !  1 < zeta, phi=5+zeta
     temp1=vkc/log(obu/z0h)
  endif
#endif

!
! Humidity profile
!

  zldis=clm%forc_hgt_q-displa
  zeta=zldis/obu
  zetat=0.465
  if (zeta < -zetat) then          ! zeta < -1
     temp2=vkc/(log(-zetat*obu/z0q) - &
          StabilityFunc(2,-zetat) + StabilityFunc(2,z0q/obu) &
          + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp2=vkc/(log(zldis/z0q) - &
          StabilityFunc(2,zeta)+StabilityFunc(2,z0q/obu))
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp2=vkc/(log(zldis/z0q)+5.*zeta-5.*z0q/obu)
  else                             !  1 < zeta, phi=5+zeta
     temp2=vkc/(log(obu/z0q) + 5. - 5.*z0q/obu &
          + (5.*log(zeta)+zeta-1.))
  endif

 !2m Calculation
  zldis=2.0 + z0q
  zeta=zldis/obu
  zetat=0.465
  if (zeta < -zetat) then          ! zeta < -1
     temp4=vkc/(log(-zetat*obu/z0q) - &
          StabilityFunc(2,-zetat) + StabilityFunc(2,z0q/obu) &
          + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp4=vkc/(log(zldis/z0q) - &
          StabilityFunc(2,zeta)+StabilityFunc(2,z0q/obu))
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp4=vkc/(log(zldis/z0q)+5.*zeta-5.*z0q/obu)
  else                             !  1 < zeta, phi=5+zeta
     temp4=vkc/(log(obu/z0q) + 5. - 5.*z0q/obu &
          + (5.*log(zeta)+zeta-1.))
  endif

#if (defined PERGRO)
  if (zeta < -zetat) then          ! zeta < -1
     temp2=vkc/log(-zetat*obu/z0q)
  else if (zeta < 0.) then         ! -1 <= zeta < 0
     temp2=vkc/log(zldis/z0q)
  else if (zeta <= 1.) then        !  0 <= ztea <= 1
     temp2=vkc/log(zldis/z0q)
  else                             !  1 < zeta, phi=5+zeta
     temp2=vkc/log(obu/z0q)
  endif
#endif

end subroutine FrictionVelocity
