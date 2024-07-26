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

subroutine MoninObukIni(ur, thv, dthv, zldis, z0m, &
                        um, obu  )

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
! Initialization of the Monin-Obukhov length.
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
! $Id: MoninObukIni.F90,v 1.6 2004/11/24 22:56:31 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2_varcon, only : grav
  implicit none

!----Arguments----------------------------------------------------------

  real(r8), intent(in) :: ur    ! wind speed at reference height [m s-1]
  real(r8), intent(in) :: thv   ! virtual potential temperature (kelvin)
  real(r8), intent(in) :: dthv  ! diff of vir. poten. temp. between ref. height and surface
  real(r8), intent(in) :: zldis ! reference height "minus" zero displacement heght [m]
  real(r8), intent(in) :: z0m   ! roughness length, momentum [m]

  real(r8), intent(out) :: um   ! wind speed including the stability effect [m s-1]
  real(r8), intent(out) :: obu  ! monin-obukhov length (m)

!----Local Variables----------------------------------------------------

  real(r8)  wc    ! convective velocity [m s-1]
  real(r8)  rib   ! bulk Richardson number
  real(r8)  zeta  ! dimensionless height used in Monin-Obukhov theory
!  real(r8)  ustar ! friction velocity [m s-1]     

!----End Variable List--------------------------------------------------

!
! Initial values of u* and convective velocity
!

!  ustar=0.06
  wc=0.5
  if (dthv >= 0.) then
     um=max(ur,0.1_r4)
  else
     um=sqrt(ur*ur+wc*wc)
  endif

  rib=grav*zldis*dthv/(thv*um*um)
#if (defined PERGRO)
  rib = 0.
#endif
  
  if (rib >= 0.) then      ! neutral or stable
     zeta = rib*log(zldis/z0m)/(1.-5.*min(rib,0.19_r4))
     zeta = min(2._r4,max(zeta,0.01_r4 ))
  else                    !unstable
     zeta=rib*log(zldis/z0m)
     zeta = max(-100._r4,min(zeta,-0.01_r4 ))
  endif

  obu=zldis/zeta

end subroutine MoninObukIni
