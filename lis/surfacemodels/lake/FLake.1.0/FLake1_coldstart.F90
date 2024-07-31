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
!BOP
!
! !ROUTINE: FLake1_coldstart
! \label{FLake1_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   6/4/13: Shugong Wang; initial implementation for LIS 7 and FLake1
!
! !INTERFACE:
subroutine FLake1_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use FLake1_Mod
   
!
! !DESCRIPTION:
!
!  This routine initializes the FLake1 state variables with some
!  predefined values uniformly for the entire domain. 
!
!EOP
 
    implicit none
    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    real    :: cb_gridDesc(6)
    do n=1, LIS_rc%nnest
        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "MSG: FLake1_coldstart -- cold-starting FLake1"
            do t=1, LIS_rc%npatch(n,mtype)
                FLAKE1_struc(n)%flake1(t)%T_snow = FLAKE1_struc(n)%init_T_snow
                FLAKE1_struc(n)%flake1(t)%T_ice = FLAKE1_struc(n)%init_T_ice
                FLAKE1_struc(n)%flake1(t)%T_mnw = FLAKE1_struc(n)%init_T_mnw
                FLAKE1_struc(n)%flake1(t)%T_wML = FLAKE1_struc(n)%init_T_wML
                FLAKE1_struc(n)%flake1(t)%T_bot = FLAKE1_struc(n)%init_T_bot
                FLAKE1_struc(n)%flake1(t)%T_b1 = FLAKE1_struc(n)%init_T_b1
                FLAKE1_struc(n)%flake1(t)%C_T = FLAKE1_struc(n)%init_C_T
                FLAKE1_struc(n)%flake1(t)%H_snow = FLAKE1_struc(n)%init_H_snow
                FLAKE1_struc(n)%flake1(t)%H_ice = FLAKE1_struc(n)%init_H_ice
                FLAKE1_struc(n)%flake1(t)%H_ML = FLAKE1_struc(n)%init_H_ML
                FLAKE1_struc(n)%flake1(t)%H_B1 = FLAKE1_struc(n)%init_H_B1
                FLAKE1_struc(n)%flake1(t)%T_sfc = FLAKE1_struc(n)%init_T_sfc
                FLAKE1_struc(n)%flake1(t)%albedo_water = FLAKE1_struc(n)%init_albedo_water
                FLAKE1_struc(n)%flake1(t)%albedo_ice = FLAKE1_struc(n)%init_albedo_ice
                FLAKE1_struc(n)%flake1(t)%albedo_snow = FLAKE1_struc(n)%init_albedo_snow
            enddo
        endif
    
        LIS_rc%yr = LIS_rc%syr
        LIS_rc%mo = LIS_rc%smo
        LIS_rc%da = LIS_rc%sda
        LIS_rc%hr = LIS_rc%shr
        LIS_rc%mn = LIS_rc%smn
        LIS_rc%ss = LIS_rc%sss
        
        call LIS_date2time(LIS_rc%time, LIS_rc%doy, LIS_rc%gmt, LIS_rc%yr,      &
                           LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
        write(LIS_logunit,*) "MSG: FLake1_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
    enddo
end subroutine FLake1_coldstart
