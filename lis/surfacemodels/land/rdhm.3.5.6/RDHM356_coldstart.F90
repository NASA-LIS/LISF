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
! !ROUTINE: RDHM356_coldstart
! \label{RDHM356_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   11/5/13: Shugong Wang; initial implementation for LIS 7 and RDHM356
!
! !INTERFACE:
subroutine RDHM356_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use RDHM356_lsmMod
   
!
! !DESCRIPTION:
!
!  This routine initializes the RDHM356 state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
 
    implicit none
    integer :: mtype
    integer :: t, l, n, i, j
    integer :: c, r
    integer :: NSOIL, NUPL, NSAC
    real :: STXT, SAND, CLAY, SUPM, SLWM, SMAX 
    real :: PSISAT, BRT, SWLT, QUARTZ, STYPE
    real :: RTUP, RTLW, ZSOIL(5), TSOIL(5)    
    real :: UZTWH, UZFWH, LZTWH, LZFSH, LZFPH
    real :: UZTWC, UZFWC, LZTWC, LZFSC, LZFPC
    real :: RSMAX, CKSL, ZBOT, TBOT, SMC(5), SH2O(5), ZX, TS0
    do n=1, LIS_rc%nnest
        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "MSG: RDHM356_coldstart -- cold-starting RDHM356"
            do t=1, LIS_rc%npatch(n,mtype)
                
                RDHM356_struc(n)%rdhm356(t)%UZTWC = RDHM356_struc(n)%init_UZTWC_ratio &
                                                  * RDHM356_struc(n)%rdhm356(t)%UZTWM
                RDHM356_struc(n)%rdhm356(t)%UZFWC = RDHM356_struc(n)%init_UZFWC_ratio &
                                                  * RDHM356_struc(n)%rdhm356(t)%UZFWM
                RDHM356_struc(n)%rdhm356(t)%LZTWC = RDHM356_struc(n)%init_LZTWC_ratio &
                                                  * RDHM356_struc(n)%rdhm356(t)%LZTWM 
                RDHM356_struc(n)%rdhm356(t)%LZFPC = RDHM356_struc(n)%init_LZFPC_ratio &
                                                  * RDHM356_struc(n)%rdhm356(t)%LZFPM 
                RDHM356_struc(n)%rdhm356(t)%LZFSC = RDHM356_struc(n)%init_LZFSC_ratio &
                                                  * RDHM356_struc(n)%rdhm356(t)%LZFSM 
                RDHM356_struc(n)%rdhm356(t)%ADIMC = RDHM356_struc(n)%init_ADIMC_ratio &
                                                  * (RDHM356_struc(n)%rdhm356(t)%UZTWM &
                                                  +  RDHM356_struc(n)%rdhm356(t)%LZTWM) 
                
                UZTWC = RDHM356_struc(n)%rdhm356(t)%UZTWC
                UZFWC = RDHM356_struc(n)%rdhm356(t)%UZFWC
                LZTWC = RDHM356_struc(n)%rdhm356(t)%LZTWC
                LZFPC = RDHM356_struc(n)%rdhm356(t)%LZFPC
                LZFSC = RDHM356_struc(n)%rdhm356(t)%LZFSC
                ! upper zone tension water storage + upper zone free water storage
                SUPM = RDHM356_struc(n)%rdhm356(t)%UZTWM + RDHM356_struc(n)%rdhm356(t)%UZFWM     
                ! lower zone tension water storage + lower zone supplemental free water storage + lower zone primary free water storage
                SLWM = RDHM356_struc(n)%rdhm356(t)%LZTWM + RDHM356_struc(n)%rdhm356(t)%LZFSM &
                     + RDHM356_struc(n)%rdhm356(t)%LZFPM 
                STXT = 1.0 * RDHM356_struc(n)%rdhm356(t)%SOILTYP 
                RSMAX = RDHM356_struc(n)%rdhm356(t)%RSMAX 
                CKSL = RDHM356_struc(n)%rdhm356(t)%CKSL
                ZBOT = RDHM356_struc(n)%rdhm356(t)%ZBOT 
                TBOT = RDHM356_struc(n)%rdhm356(t)%TBOT 
                SAND = RDHM356_struc(n)%rdhm356(t)%SAND 
                CLAY = RDHM356_struc(n)%rdhm356(t)%CLAY 

                call SOILPAR2(STXT, SAND, CLAY, SUPM, SLWM, SMAX,   & 
                              PSISAT, BRT, SWLT, QUARTZ, STYPE,     &
                              NSOIL, NUPL, NSAC, ZSOIL, RTUP, RTLW)    
                
                ! intialize soil temperature according to ts0 and tbot 
                TS0 = RDHM356_struc(n)%init_TS0
                do j=1, NSOIL
                    if(j .eq. 1) then
                        ZX = -0.5 * ZSOIL(1)
                    else
                        ZX = -0.5 * (ZSOIL(j-1)+ZSOIL(j))
                    endif
                    TSOIL(j) = TS0 + (TBOT-273.16-TS0)*ZX/(-1.0*ZBOT)
                enddo

                RDHM356_struc(n)%rdhm356(t)%TS0 = TSOIL(1) 
                RDHM356_struc(n)%rdhm356(t)%TS1 = TSOIL(2)
                RDHM356_struc(n)%rdhm356(t)%TS2 = TSOIL(3)
                RDHM356_struc(n)%rdhm356(t)%TS3 = TSOIL(4)
                RDHM356_struc(n)%rdhm356(t)%TS4 = TSOIL(5)

                call INIT_SOIL_MOIST(NSOIL, NUPL,                                   &
                                     UZTWC, UZFWC, LZTWC, LZFSC, LZFPC,             &
                                     RSMAX, CKSL, ZBOT, RTUP, RTLW, PSISAT, SWLT,   &
                                     TSOIL, ZSOIL, TBOT, BRT, SMAX,                 &
                                     UZTWH, UZFWH, LZTWH, LZFSH, LZFPH, SMC, SH2O)

                RDHM356_struc(n)%rdhm356(t)%UZTWH = UZTWH
                RDHM356_struc(n)%rdhm356(t)%UZFWH = UZFWH
                RDHM356_struc(n)%rdhm356(t)%LZTWH = LZTWH
                RDHM356_struc(n)%rdhm356(t)%LZFSH = LZFSH
                RDHM356_struc(n)%rdhm356(t)%LZFPH = LZFPH

                do l=1, 5 ! TODO: check loop
                    RDHM356_struc(n)%rdhm356(t)%SMC(l) = SMC(l)
                enddo
                RDHM356_struc(n)%rdhm356(t)%SMC(6) = 0.0 
                do l=1, 5 ! TODO: check loop
                    RDHM356_struc(n)%rdhm356(t)%SH2O(l) = SH2O(l)
                enddo
                RDHM356_struc(n)%rdhm356(t)%SH2O(6) = 0.0 

                RDHM356_struc(n)%rdhm356(t)%WE = RDHM356_struc(n)%init_WE
                RDHM356_struc(n)%rdhm356(t)%LIQW = RDHM356_struc(n)%init_LIQW
                RDHM356_struc(n)%rdhm356(t)%NEGHS = RDHM356_struc(n)%init_NEGHS
                RDHM356_struc(n)%rdhm356(t)%TINDEX = RDHM356_struc(n)%init_TINDEX
                RDHM356_struc(n)%rdhm356(t)%ACCMAX = RDHM356_struc(n)%init_ACCMAX
                RDHM356_struc(n)%rdhm356(t)%SNDPT = RDHM356_struc(n)%init_SNDPT
                RDHM356_struc(n)%rdhm356(t)%SNTMP = RDHM356_struc(n)%init_SNTMP
                RDHM356_struc(n)%rdhm356(t)%SB = RDHM356_struc(n)%init_SB
                RDHM356_struc(n)%rdhm356(t)%SBAESC = RDHM356_struc(n)%init_SBAESC
                RDHM356_struc(n)%rdhm356(t)%SBWS = RDHM356_struc(n)%init_SBWS
                RDHM356_struc(n)%rdhm356(t)%STORAGE = RDHM356_struc(n)%init_STORAGE
                RDHM356_struc(n)%rdhm356(t)%AEADJ = RDHM356_struc(n)%init_AEADJ
                do l=1, 7 ! TODO: check loop
                    RDHM356_struc(n)%rdhm356(t)%EXLAG(l) = RDHM356_struc(n)%init_EXLAG(l)
                enddo
                RDHM356_struc(n)%rdhm356(t)%NEXLAG = RDHM356_struc(n)%init_NEXLAG
                RDHM356_struc(n)%rdhm356(t)%TA_PREV = RDHM356_struc(n)%init_TA_PREV
                RDHM356_struc(n)%rdhm356(t)%CH = 1e-4;
                RDHM356_struc(n)%rdhm356(t)%CM = 1e-4; 
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
        write(LIS_logunit,*) "MSG: RDHM356_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
    enddo
end subroutine RDHM356_coldstart
