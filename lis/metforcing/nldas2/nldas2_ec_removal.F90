!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: nldas2_ec_removal
! \label{nldas2_ec_removal}
!
! !REVISION HISTORY:
!  20 Oct 2006: Kristi Arsenault; Adapted elevation correction code
!               to remove such a correction from NLDAS2 fields
!
! !INTERFACE:

subroutine nldas2_ec_removal( nest, point, force_tmp, force_hum, &
                              force_lwd, force_prs )

! !USES:
  use LIS_coreMod,        only : LIS_rc
  use nldas2_forcingMod,  only : nldas2_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: nest

! !DESCRIPTION:
!  Removes Temperature, Pressure, Humidity and Longwave Radiation
!  forcing correction.
!
!  The corrections are based on the lapse-rate and hypsometric adjustments 
!  to these variables described in Cosgrove et. al (2003).
!  
!  Cosgrove, B.A. et.al, Real-time and retrospective forcing in the 
!  North American Land Data Assimilation (NLDAS2) project, Journal of 
!  Geophysical Research, 108(D22), 8842, DOI: 10.1029/2002JD003118, 2003.  
!
!  The arguments are: 
!  \begin{description}
!   \item [nest]
!     index of the domain or nest.
!   \item [point]
!     index of the grid point
!   \item [force\_tmp]
!     temperature value for the grid point
!   \item [force\_hum]
!     specific humidity for the grid point
!   \item [force\_lwd]
!     downward longwave radiation for the grid point
!   \item [force\_prs]
!    surface pressure for the grid point
!  \end{description}
!
!EOP

  integer, intent (in) :: point
  real, intent (inout) :: force_tmp, force_hum
  real, intent (inout) :: force_lwd, force_prs
  real :: orig_tmp, orig_hum
  real :: orig_lwd, orig_prs
  real :: elevdiff

  real :: mee, mfe, ee, fe, ratio
  real :: esat,qsat,rh_corr,fesat,fqsat,femiss,emiss
  real :: tbar

  integer, parameter :: bb = 2016
  real, parameter    :: grav = 9.81
  real, parameter    :: rdry = 287.
  real, parameter    :: lapse = -0.0065

! ----------------------------------------------------------------

     elevdiff = nldas2_struc(nest)%orig_ediff(point)

! -- Apply elevation correction to temperature::
     orig_tmp = force_tmp - (lapse*elevdiff)
     tbar = (orig_tmp + force_tmp) / 2.

! -- Apply elevation correction to surface pressure::
     orig_prs = force_prs*(exp((grav*elevdiff)/(rdry*tbar)))

! -- Apply elevation correction to humidity::

     if (force_hum .eq. 0) force_hum=1e-08

     esat = 611.2*(exp((17.67*(orig_tmp-273.15))/ &
                       ((orig_tmp-273.15)+243.5)))
     qsat = (0.622*esat)/(orig_prs-(0.378*esat))

     fesat=611.2*(exp((17.67*(force_tmp-273.15))/ &
                     ((force_tmp-273.15)+243.5)))

     fqsat=(0.622*fesat)/(force_prs-(0.378*fesat))
     rh_corr = (force_hum / fqsat) * 100.

     orig_hum = (rh_corr * qsat) / 100.

! -- Apply elevation correction to downward LW radiation::
     ee = (orig_hum*orig_prs) / 0.622
     fe = (force_hum*force_prs) / 0.622
     mee = ee/100.
     mfe = fe/100.

!----------------------------------------------------------------------
! Correct for negative vapor pressure at very low temperatures at
!  high latitudes
!----------------------------------------------------------------------
     if (mee .le. 0) mee = 1e-08
     if (mfe .le. 0) mfe = 1e-08

     emiss = 1.08*(1-exp(-mee**(orig_tmp/bb)))  
     femiss = 1.08*(1-exp(-mfe**(force_tmp/bb)))
     ratio = (femiss*(force_tmp**4))/(emiss*(orig_tmp**4))

     orig_lwd = force_lwd / ratio

!-- Reassign uncorrected fields to LIS forcing fields::
    force_tmp = orig_tmp
    force_hum = orig_hum
    force_lwd = orig_lwd
    force_prs = orig_prs


end subroutine nldas2_ec_removal
