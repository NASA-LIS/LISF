!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: Ac71_coldstart
! \label{Ac71_coldstart}
!
! !REVISION HISTORY: 
!   18 JAN 2024, Louise Busschaert; initial implementation for LIS 7 and AC71
!
! !INTERFACE:
subroutine Ac71_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use Ac71_lsmMod
 !  !AC71
    use ac_global, only:    GetCompartment_Layer, &
                            GetCompartment_theta         
!
! !DESCRIPTION:
!
!  This routine initializes the AC71 state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
 
    implicit none
    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    
    do n=1, LIS_rc%nnest
        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "MSG: Ac71_coldstart -- cold-starting Ac71"
            do t=1, LIS_rc%npatch(n,mtype)
                !AC71_struc(n)%ac71(t)%irun = 1
                !AC71_struc(n)%daynrinextclimaterecord = 1

                !AC71_struc(n)%ac71(t)%InitializeRun = 1 ! gets 1 at end of year 

                do l=1, AC71_struc(n)%ac71(t)%NrCompartments
                    AC71_struc(n)%ac71(t)%smc(l) = AC71_struc(n)%init_smc(GetCompartment_Layer(l))
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
        write(LIS_logunit,*) "MSG: Ac71_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
    enddo
end subroutine Ac71_coldstart
