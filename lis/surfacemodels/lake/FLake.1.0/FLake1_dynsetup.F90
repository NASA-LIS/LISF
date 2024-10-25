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
! !ROUTINE: FLake1_dynsetup
! \label{FLake1_dynsetup}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   6/4/13: Shugong Wang; initial implementation for LIS 7 and FLake1
!
! !INTERFACE:
subroutine FLake1_dynsetup(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_domain, LIS_surface
    use LIS_timeMgrMod, only : LIS_date2time, LIS_tick
    use FLake1_Mod, only : FLAKE1_struc
   !use any other modules 
!
! !DESCRIPTION:
!  This routine sets up the time-dependent variables in FLake1
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
    integer   :: col, row, ncount(LIS_rc%npatch(n, LIS_rc%lake_index))
!TODO: add code here
 
end subroutine FLake1_dynsetup
