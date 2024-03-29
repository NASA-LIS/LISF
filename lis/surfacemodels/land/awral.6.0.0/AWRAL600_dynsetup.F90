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
! !ROUTINE: AWRAL600_dynsetup
! \label{AWRAL600_dynsetup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   12/18/18: Wendy Sharples, Shugong Wang; initial implementation for LIS 7 and AWRAL600
!
! !INTERFACE:
subroutine AWRAL600_dynsetup(n)
! !USES:
    use LIS_logMod, only     : LIS_logunit
    use LIS_coreMod, only    : LIS_rc, LIS_domain, LIS_surface
    use LIS_timeMgrMod, only : LIS_date2time, LIS_tick
    use LIS_fileIOMod, only  : LIS_read_param
    use AWRAL600_lsmMod, only : AWRAL600_struc
   !use any other modules 
!
! !DESCRIPTION:
!  This routine sets up the time-dependent variables in AWRAL600
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
end subroutine AWRAL600_dynsetup
 
! generate date/time string for reading time-dependent variables 
