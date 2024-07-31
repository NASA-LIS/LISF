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

subroutine Hydrology_Lake (clm) 

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
!
! Method:
!
! Author:
! Gordon Bonan
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!
!-----------------------------------------------------------------------
! $Id: Hydrology_Lake.F90,v 1.6 2004/11/24 22:56:28 jim Exp $
!-----------------------------------------------------------------------

  use LIS_precisionMod
  use clm2type
  use clm2_varpar, only : nlevsoi
  use clm2_varcon, only : hfus, tfrz, spval, denh2o
  implicit none

!----Arguments----------------------------------------------------------

  type (clm1d), intent(inout) :: clm      ! CLM 1-D Module

!----Local Variables----------------------------------------------------

  integer j               ! do loop index
  real(r8) qflx_evap_grnd ! ground surface evaporation rate (mm h2o/s)
  real(r8) qflx_dew_grnd  ! ground surface dew formation (mm h2o /s) [+]
  real(r8) qflx_sub_snow  ! sublimation rate from snow pack (mm h2o /s) [+]
  real(r8) qflx_dew_snow  ! surface dew added to snow pack (mm h2o /s) [+]

!----End Variable List--------------------------------------------------

!
! [1] Snow on the lake ice 
!
  qflx_evap_grnd = 0.
  qflx_sub_snow = 0.
  qflx_dew_snow = 0.
  qflx_dew_grnd = 0.

  if (clm%qflx_evap_soi >= 0.) then

!
! Sublimation: do not allow for more sublimation than there is snow
! after melt.  Remaining surface evaporation used for infiltration.
!

     qflx_sub_snow = min( clm%qflx_evap_soi, clm%h2osno/clm%dtime-clm%qmelt )
     qflx_evap_grnd = clm%qflx_evap_soi - qflx_sub_snow

  else

     if (clm%t_grnd < tfrz-0.1) then
        qflx_dew_snow = abs(clm%qflx_evap_soi)
     else
        qflx_dew_grnd = abs(clm%qflx_evap_soi)
     endif

  endif

!
! Update snow pack
!

  if (clm%do_capsnow) then 
     clm%h2osno = clm%h2osno - (clm%qmelt + qflx_sub_snow)*clm%dtime
  else
     clm%h2osno = clm%h2osno + (clm%forc_snow-clm%qmelt-qflx_sub_snow+qflx_dew_snow)*clm%dtime
  endif
  clm%h2osno = max( clm%h2osno, 0._r4 )

!
! No snow if lake unfrozen
!

  if (clm%t_grnd > tfrz) clm%h2osno = 0.

!
! Snow depth
!

  clm%snowdp = clm%h2osno/250.    !Assume a constant snow bulk density = 250.

!
! Determine ending water balance
!
  
  clm%endwb = clm%h2osno

!
! [2] Set other clm values for lake points
!

!
! The following are needed for global average on history tape. 
! Note: time invariant variables set in initialization phase:
! z, dz, snl, h2osoi_liq, and h2osoi_ice 
!

  clm%eflx_snomelt    = clm%qmelt*hfus
  clm%h2ocan          = 0.  
  clm%qflx_evap_veg   = 0.
  clm%qflx_tran_veg   = 0.
  clm%qflx_infl       = 0. 
  clm%qflx_snomelt    = clm%qmelt
  clm%qflx_surf       = 0.
  clm%qflx_drain      = 0.
  clm%qflx_qrgwl      = clm%forc_rain + clm%forc_snow - clm%qflx_evap_tot - &
                        (clm%endwb-clm%begwb)/clm%dtime
  clm%qflx_prec_grnd  = clm%forc_rain + clm%forc_snow
  clm%qflx_prec_intr  = 0.

! Components that are not displayed over lake on history tape and 
! therefore need to be set to spval here
#if 0 
  clm%btran         = spval
  clm%rootr(:)      = spval
  clm%snowice       = spval
  clm%snowliq       = spval
  clm%h2osoi_vol(:) = spval
  clm%h2osoi_ice(:) = spval
  clm%h2osoi_liq(:) = spval
#endif
  clm%btran         = spval
  clm%rootfr(:)      = 1.0
  clm%snowice       = 0
  clm%snowliq       = 0
  clm%h2osoi_vol(:) = 0.5
  clm%h2osoi_ice(:) = 0.0
  do j=1,nlevsoi
     clm%h2osoi_liq(j) = 0.5*clm%dz(j)*denh2o
  enddo   
end subroutine Hydrology_Lake
