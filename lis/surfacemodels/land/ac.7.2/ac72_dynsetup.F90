!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: AC72_dynsetup
! \label{AC72_dynsetup}
!
! !REVISION HISTORY:
!  18 JAN 2024, Louise Busschaert; initial implementation for AC72
!
! !INTERFACE:
subroutine AC72_dynsetup(n)
! !USES:
    use LIS_logMod, only     : LIS_logunit
    use LIS_coreMod, only    : LIS_rc, LIS_domain, LIS_surface
    use LIS_timeMgrMod, only : LIS_date2time, LIS_tick
    use LIS_fileIOMod, only  : LIS_read_param
    use AC72_lsmMod, only : AC72_struc
!
! !DESCRIPTION:
!  This routine sets up the time-dependent variables in AC72
!
!EOP
    implicit none
    integer, intent(in) :: n
 
    !TODO: add code here if needed.
end subroutine AC72_dynsetup
 
! generate date/time string for reading time-dependent variables ??Michel
