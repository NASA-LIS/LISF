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
! !ROUTINE: AWRAL600_coldstart
! \label{AWRAL600_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   12/18/18: Wendy Sharples, Shugong Wang; initial implementation for LIS 7 and AWRAL600
!
! !INTERFACE:
subroutine AWRAL600_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use AWRAL600_lsmMod
   
!
! !DESCRIPTION:
!
!  This routine initializes the AWRAL600 state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
 
    implicit none
    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    do n=1, LIS_rc%nnest
        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "[INFO]: AWRAL600_coldstart -- cold-starting AWRAL600"
            do t=1, LIS_rc%npatch(n,mtype)
                AWRAL600_struc(n)%awral600(t)%sr = AWRAL600_struc(n)%init_sr
                AWRAL600_struc(n)%awral600(t)%sg = AWRAL600_struc(n)%init_sg
                do l=1, AWRAL600_struc(n)%nhru
                    AWRAL600_struc(n)%awral600(t)%s0(l) = AWRAL600_struc(n)%init_s0(l)
                enddo
                do l=1, AWRAL600_struc(n)%nhru
                    AWRAL600_struc(n)%awral600(t)%ss(l) = AWRAL600_struc(n)%init_ss(l)
                enddo
                do l=1, AWRAL600_struc(n)%nhru
                    AWRAL600_struc(n)%awral600(t)%sd(l) = AWRAL600_struc(n)%init_sd(l)
                enddo
                do l=1, AWRAL600_struc(n)%nhru
                    AWRAL600_struc(n)%awral600(t)%mleaf(l) = AWRAL600_struc(n)%init_mleaf(l)
                enddo
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
        write(LIS_logunit,*) "[INFO]: AWRAL600_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
    enddo
end subroutine AWRAL600_coldstart
