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
! !ROUTINE: RDHM356_readrst
! \label{RDHM356_readrst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   3/6/14: Shugong Wang; initial implementation for LIS 7 and RDHM356
!
! !INTERFACE:
subroutine RDHM356_readrst()
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify                
    use RDHM356_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for RDHM356.  This
!  includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the RDHM356
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    UZTWC                      - RDHM356 upper zone tension water storage content [mm]
!    UZFWC                      - RDHM356 upper zone free water storage content [mm]
!    LZTWC                      - RDHM356 lower zone tension water storage content [mm]
!    LZFPC                      - RDHM356 lower zone primary free water storage content [mm]
!    LZFSC                      - RDHM356 lower zone supplemental free water storage content [mm]
!    ADIMC                      - RDHM356 additional impervious area content [mm]
!    TS0                        - RDHM356 first soil layer temperature [C]
!    TS1                        - RDHM356 second soil layer temperature [C]
!    TS2                        - RDHM356 third soil layer temperature [C]
!    TS3                        - RDHM356 fourth soil layer temperature [C]
!    TS4                        - RDHM356 fifth soil layer temperature [C]
!    UZTWH                      - RDHM356 unfrozen upper zone tension water [mm]
!    UZFWH                      - RDHM356 unfrozen uppeer zone free water [mm]
!    LZTWH                      - RDHM356 unfrozen lower zone tension water [mm]
!    LZFSH                      - RDHM356 unfrozen lower zone supplemental free water [mm]
!    LZFPH                      - RDHM356 unfrozen lower zone primary free water [mm]
!    SMC                        - RDHM356 volumetric content of total soil moisture at each layer [m^3 m-3]
!    SH2O                       - RDHM356 volumetric content of liquid soil moisture at each layer [m^3 m-3]
!    WE                         - RDHM356 snow water equivalent without liquid water [mm]
!    LIQW                       - RDHM356 liquid water in snow [mm]
!    NEGHS                      - RDHM356 negative snow heat [mm]
!    TINDEX                     - RDHM356 antecedent temperature index [C]
!    ACCMAX                     - RDHM356 cumulated snow water including liquid [mm]
!    SNDPT                      - RDHM356 snow depth [cm]
!    SNTMP                      - RDHM356 average snow temperature [C]
!    SB                         - RDHM356 the last highest snow water equivalent before any snow fall [C]
!    SBAESC                     - RDHM356 internal snow state during melt & new snow fall (checked with Victor) [-]
!    SBWS                       - RDHM356 internal snow state during melt & new snow fall (checked with Victor) [-]
!    STORAGE                    - RDHM356 snow liquid water attenuation storage [mm]
!    AEADJ                      - RDHM356 adjusted areal snow cover fraction [-]
!    EXLAG                      - RDHM356 array of lagged liquid water values [-]
!    NEXLAG                     - RDHM356 number of ordinates in lagged liquid water array (EXLAG) [-]
!    TA_PREV                    - RDHM356 air temperature of previous time step [-]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!      reads a variable from the restart file
!   \item[RDHM356\_coldstart](\ref{RDHM356_coldstart}) \newline
!      initializes the RDHM356 state variables
! \end{description}
!EOP
 
    implicit none
 
    integer           :: t, l
    integer           :: nc, nr, npatch
    integer           :: n
    integer           :: ftn
    integer           :: status
    real, allocatable :: tmptilen(:)
    logical           :: file_exists
    character*20      :: wformat
 
    do n=1, LIS_rc%nnest
        wformat = trim(RDHM356_struc(n)%rformat)
        ! coldstart
        if(LIS_rc%startcode .eq. "coldstart") then  
            call RDHM356_coldstart(LIS_rc%lsm_index)
        ! restart
        elseif(LIS_rc%startcode .eq. "restart") then
            allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            ! check the existance of restart file
            inquire(file=RDHM356_struc(n)%rfile, exist=file_exists)
            If (.not. file_exists) then 
                write(LIS_logunit,*) "RDHM356 restart file ", RDHM356_struc(n)%rfile," does not exist "
                write(LIS_logunit,*) "Program stopping ..."
                call LIS_endrun
            endif
            write(LIS_logunit,*) "RDHM356 restart file used: ", RDHM356_struc(n)%rfile
        
            ! open restart file
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn, file=RDHM356_struc(n)%rfile, form="unformatted")
                read(ftn) nc, nr, npatch  !time, veg class, no. tiles
 
                ! check for grid space conflict
                if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
                    write(LIS_logunit,*) RDHM356_struc(n)%rfile, "grid space mismatch - RDHM356 halted"
                    call LIS_endrun
                endif
            
                if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
                    write(LIS_logunit,*) "restart tile space mismatch, halting..."
                    call LIS_endrun
                endif
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_open(path=RDHM356_struc(n)%rfile, &
                                   mode=NF90_NOWRITE, ncid=ftn)
                call LIS_verify(status, "Error opening file "//RDHM356_struc(n)%rfile)
#endif
            endif
 
            ! read: upper zone tension water storage content
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%UZTWC, &
                                     varname="UZTWC", wformat=wformat)
 
            ! read: upper zone free water storage content
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%UZFWC, &
                                     varname="UZFWC", wformat=wformat)
 
            ! read: lower zone tension water storage content
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZTWC, &
                                     varname="LZTWC", wformat=wformat)
 
            ! read: lower zone primary free water storage content
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZFPC, &
                                     varname="LZFPC", wformat=wformat)
 
            ! read: lower zone supplemental free water storage content
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZFSC, &
                                     varname="LZFSC", wformat=wformat)
 
            ! read: additional impervious area content
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%ADIMC, &
                                     varname="ADIMC", wformat=wformat)
 
            ! read: first soil layer temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TS0, &
                                     varname="TS0", wformat=wformat)
 
            ! read: second soil layer temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TS1, &
                                     varname="TS1", wformat=wformat)
 
            ! read: third soil layer temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TS2, &
                                     varname="TS2", wformat=wformat)
 
            ! read: fourth soil layer temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TS3, &
                                     varname="TS3", wformat=wformat)
 
            ! read: fifth soil layer temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TS4, &
                                     varname="TS4", wformat=wformat)
 
            ! read: unfrozen upper zone tension water
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%UZTWH, &
                                     varname="UZTWH", wformat=wformat)
 
            ! read: unfrozen uppeer zone free water
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%UZFWH, &
                                     varname="UZFWH", wformat=wformat)
 
            ! read: unfrozen lower zone tension water
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZTWH, &
                                     varname="LZTWH", wformat=wformat)
 
            ! read: unfrozen lower zone supplemental free water
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZFSH, &
                                     varname="LZFSH", wformat=wformat)
 
            ! read: unfrozen lower zone primary free water
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZFPH, &
                                     varname="LZFPH", wformat=wformat)
 
            ! read: volumetric content of total soil moisture at each layer
            do l=1, 6 ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SMC", &
                                         dim=l, vlevels = 6, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    RDHM356_struc(n)%rdhm356(t)%SMC(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: volumetric content of liquid soil moisture at each layer
            do l=1, 6 ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SH2O", &
                                         dim=l, vlevels = 6, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    RDHM356_struc(n)%rdhm356(t)%SH2O(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: snow water equivalent without liquid water
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%WE, &
                                     varname="WE", wformat=wformat)
 
            ! read: liquid water in snow
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LIQW, &
                                     varname="LIQW", wformat=wformat)
 
            ! read: negative snow heat
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%NEGHS, &
                                     varname="NEGHS", wformat=wformat)
 
            ! read: antecedent temperature index
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TINDEX, &
                                     varname="TINDEX", wformat=wformat)
 
            ! read: cumulated snow water including liquid
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%ACCMAX, &
                                     varname="ACCMAX", wformat=wformat)
 
            ! read: snow depth
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%SNDPT, &
                                     varname="SNDPT", wformat=wformat)
 
            ! read: average snow temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%SNTMP, &
                                     varname="SNTMP", wformat=wformat)
 
            ! read: the last highest snow water equivalent before any snow fall
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%SB, &
                                     varname="SB", wformat=wformat)
 
            ! read: internal snow state during melt & new snow fall (checked with Victor)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%SBAESC, &
                                     varname="SBAESC", wformat=wformat)
 
            ! read: internal snow state during melt & new snow fall (checked with Victor)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%SBWS, &
                                     varname="SBWS", wformat=wformat)
 
            ! read: snow liquid water attenuation storage
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%STORAGE, &
                                     varname="STORAGE", wformat=wformat)
 
            ! read: adjusted areal snow cover fraction
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%AEADJ, &
                                     varname="AEADJ", wformat=wformat)
 
            ! read: array of lagged liquid water values
            do l=1, 7 ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="EXLAG", &
                                         dim=l, vlevels = 7, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    RDHM356_struc(n)%rdhm356(t)%EXLAG(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: number of ordinates in lagged liquid water array (EXLAG)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%NEXLAG, &
                                     varname="NEXLAG", wformat=wformat)
 
            ! read: air temperature of previous time step
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TA_PREV, &
                                     varname="TA_PREV", wformat=wformat)
            ! read: ch 
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%CH, &
                                     varname="CH", wformat=wformat)
        
            ! read: cm 
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%CM, &
                                     varname="CM", wformat=wformat)
            ! close restart file
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in RDHM356_readrst")
#endif
            endif
            deallocate(tmptilen)
        endif    
    enddo
end subroutine RDHM356_readrst
        
