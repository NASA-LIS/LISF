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

subroutine SensibleHCond (ra,   rb,   rd,   wta,  wtl,  &
                          wtg,  wta0, wtl0, wtg0, wtal, &
                          wtga, wtgl, clm   ) 

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
! Provides dimensional and non-dimensional sensible heat
! conductances for canopy and soil flux calculations.
!
! Method:
!
! Author:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: SensibleHCond.F90,v 1.6 2004/11/24 22:56:34 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout)  :: clm !CLM 1-D Module

  real(r8), intent(in) :: ra ! aerodynamical resistance [s/m]
  real(r8), intent(in) :: rb ! leaf boundary layer resistance [s/m]
  real(r8), intent(in) :: rd ! thermal resistance between ground and bottom of canopy

  real(r8), intent(out) :: wta  ! heat conduactance for air [m s-1]
  real(r8), intent(out) :: wtg  ! heat conduactance for ground [m s-1]
  real(r8), intent(out) :: wtl  ! heat conduactance for leaf [m s-1]
  real(r8), intent(out) :: wta0 ! normalized heat conductance for air [-]
  real(r8), intent(out) :: wtl0 ! normalized heat conductance for air [-]
  real(r8), intent(out) :: wtg0 ! normalized heat conductance for ground [-]
  real(r8), intent(out) :: wtal ! normalized heat conductance for air and leaf [-]
  real(r8), intent(out) :: wtgl ! normalized heat conductance for leaf and ground [-]
  real(r8), intent(out) :: wtga ! normalized heat conductance for air and ground  [-]

!----Local Variables----------------------------------------------------

  real(r8) wtshi                ! heat resistance for air, ground and leaf [s/m]

!----End Variable List--------------------------------------------------

  wta   = 1./ra                     ! air
  wtl   = (clm%elai+clm%esai)/rb    ! leaf
  wtg   = 1./rd                     ! ground
  wtshi = 1./(wta+wtl+wtg)

  wtl0  = wtl*wtshi         ! leaf
  wtg0  = wtg*wtshi         ! ground
  wta0  = wta*wtshi         ! air

  wtgl  = wtl0+wtg0         ! ground + leaf
  wtga  = wta0+wtg0         ! ground + air
  wtal  = wta0+wtl0         ! air + leaf

end subroutine SensibleHCond
