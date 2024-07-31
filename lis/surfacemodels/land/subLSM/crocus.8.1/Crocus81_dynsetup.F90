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
!
! !ROUTINE: Crocus81_dynsetup
! \label{Crocus81_dynsetup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/18/19: Mahdi Navari, Shugong Wang; initial implementation for LIS 7 and Crocus81
!
! !INTERFACE:
subroutine Crocus81_dynsetup(n)
! !USES:
    use LIS_logMod, only     : LIS_logunit
    use LIS_coreMod, only    : LIS_rc, LIS_domain, LIS_surface
    use LIS_timeMgrMod, only : LIS_date2time, LIS_tick
    use LIS_fileIOMod, only  : LIS_read_param
    use Crocus81_lsmMod, only : CROCUS81_struc
   !use any other modules 
!
! !DESCRIPTION:
!  This routine sets up the time-dependent variables in Crocus81
!
!EOP
    implicit none
    integer, intent(in) :: n
    
    integer   :: tid
    integer   :: t, gid, change, local_hour
    integer   :: locdoy, locyr, locmo, locda, lochr, locmn, locss
    real*8    :: loctime
    real      :: interp_fraction
    real      :: locgmt
    integer   :: col, row, ncount(LIS_rc%npatch(n, LIS_rc%lsm_index))
 
    !TODO: add code here if needed.
end subroutine Crocus81_dynsetup
 
! generate date/time string for reading time-dependent variables 
