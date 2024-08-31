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
! !ROUTINE: noahmpglacier3911_dynsetup
! \label{noahmpglacier3911_dynsetup}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine noahmpglacier3911_dynsetup(n)
! !USES:
    use LIS_logMod, only     : LIS_logunit
    use LIS_coreMod, only    : LIS_rc, LIS_domain, LIS_surface
    use LIS_timeMgrMod, only : LIS_date2time, LIS_tick
    use LIS_fileIOMod, only  : LIS_read_param
   !use any other modules 
!
! !DESCRIPTION:
!  This routine sets up the time-dependent variables in noahmpglacier3911
!
!EOP
    implicit none
    integer, intent(in) :: n

end subroutine noahmpglacier3911_dynsetup
 
! generate date/time string for reading time-dependent variables 
