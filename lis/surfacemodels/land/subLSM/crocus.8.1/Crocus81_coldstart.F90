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
! !ROUTINE: Crocus81_coldstart
! \label{Crocus81_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/18/19: Mahdi Navari, Shugong Wang; initial implementation for LIS 7 and Crocus81
!
! !INTERFACE:
subroutine Crocus81_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use Crocus81_lsmMod
   
!
! !DESCRIPTION:
!
!  This routine initializes the Crocus81 state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
 
    implicit none
    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    do n=1, LIS_rc%nnest
        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "MSG: Crocus81_coldstart -- cold-starting Crocus81"
            do t=1, LIS_rc%npatch(n,mtype)
                do l=1, CROCUS81_struc(n)%nsnow
                    CROCUS81_struc(n)%crocus81(t)%SNOWSWE(l) = CROCUS81_struc(n)%init_SNOWSWE(l)
                enddo
                do l=1, CROCUS81_struc(n)%nsnow
                    CROCUS81_struc(n)%crocus81(t)%SNOWRHO(l) = CROCUS81_struc(n)%init_SNOWRHO(l)
                enddo
                do l=1, CROCUS81_struc(n)%nsnow
                    CROCUS81_struc(n)%crocus81(t)%SNOWHEAT(l) = CROCUS81_struc(n)%init_SNOWHEAT(l)
                enddo
                CROCUS81_struc(n)%crocus81(t)%SNOWALB = CROCUS81_struc(n)%init_SNOWALB
                do l=1, CROCUS81_struc(n)%nsnow
                    CROCUS81_struc(n)%crocus81(t)%SNOWGRAN1(l) = CROCUS81_struc(n)%init_SNOWGRAN1(l)
                enddo
                do l=1, CROCUS81_struc(n)%nsnow
                    CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2(l) = CROCUS81_struc(n)%init_SNOWGRAN2(l)
                enddo
                do l=1, CROCUS81_struc(n)%nsnow
                    CROCUS81_struc(n)%crocus81(t)%SNOWHIST(l) = CROCUS81_struc(n)%init_SNOWHIST(l)
                enddo
                do l=1, CROCUS81_struc(n)%nsnow
                    CROCUS81_struc(n)%crocus81(t)%SNOWAGE(l) = CROCUS81_struc(n)%init_SNOWAGE(l)
                enddo
                do l=1, CROCUS81_struc(n)%nsnow
                    CROCUS81_struc(n)%crocus81(t)%SNOWLIQ(l) = CROCUS81_struc(n)%init_SNOWLIQ(l)
                enddo
                do l=1, CROCUS81_struc(n)%nsnow
                    CROCUS81_struc(n)%crocus81(t)%SNOWTEMP(l) = CROCUS81_struc(n)%init_SNOWTEMP(l)
                enddo
                do l=1, CROCUS81_struc(n)%nsnow
                    CROCUS81_struc(n)%crocus81(t)%SNOWDZ(l) = CROCUS81_struc(n)%init_SNOWDZ(l)
                enddo
                CROCUS81_struc(n)%crocus81(t)%GRNDFLUX = CROCUS81_struc(n)%init_GRNDFLUX
                CROCUS81_struc(n)%crocus81(t)%SNDRIFT = CROCUS81_struc(n)%init_SNDRIFT
                CROCUS81_struc(n)%crocus81(t)%RI_n = CROCUS81_struc(n)%init_RI_n
                CROCUS81_struc(n)%crocus81(t)%CDSNOW = CROCUS81_struc(n)%init_CDSNOW
                CROCUS81_struc(n)%crocus81(t)%USTARSNOW = CROCUS81_struc(n)%init_USTARSNOW
                CROCUS81_struc(n)%crocus81(t)%CHSNOW = CROCUS81_struc(n)%init_CHSNOW
                CROCUS81_struc(n)%crocus81(t)%SNOWMAK_dz = CROCUS81_struc(n)%init_SNOWMAK_dz
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
        write(LIS_logunit,*) "[INFO]: Crocus81_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
    enddo
end subroutine Crocus81_coldstart
