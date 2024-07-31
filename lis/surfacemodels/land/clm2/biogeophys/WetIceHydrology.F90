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

subroutine WetIceHydrology (clm)

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
! Calculate hydrology for ice and wetland
!
! Method:
!
! Author:
! 7 November 2000: Mariana Vertenstein; Initial code
!
!-----------------------------------------------------------------------
! $Id: WetIceHydrology.F90,v 1.6 2004/11/24 22:56:50 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varcon, only : istwet, istice
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm !CLM 1-D Module

!----Local Variables----------------------------------------------------
! None
!----End Variable List--------------------------------------------------

!
! Wetland and land ice runoff
!

  clm%qflx_drain  = 0. 
  clm%qflx_surf   = 0.
  clm%qflx_infl   = 0.
  clm%qflx_qrgwl  = clm%forc_rain + clm%forc_snow - clm%qflx_evap_tot - &
                    (clm%endwb - clm%begwb)/clm%dtime

end subroutine WetIceHydrology
