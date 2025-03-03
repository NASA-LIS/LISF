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
! !ROUTINE: NoahMP36_readrst
! \label{NoahMP36_readrst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial
!  specification of the subroutine is defined by Sujay Kumar.
!   9/4/14: Shugong Wang; initial implementation for LIS 7 and NoahMP36
!
! !INTERFACE:
subroutine NoahMP36_readrst()
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify
    use NoahMP36_lsmMod
    !WN
    use ESMF
    use LIS_fileIOMod
    use LIS_timeMgrMod
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    !-----------------
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for NoahMP36.  This
!  includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the NoahMP36
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    albold                     - NoahMP36 snow albedo at last time step [-]
!    sneqvo                     - NoahMP36 snow mass at the last time step [mm]
!    sstc                       - NoahMP36 snow/soil temperature [K]
!    sh2o                       - NoahMP36 volumetric liquid soil moisture [m^3 m-3]
!    smc                        - NoahMP36 volumetric soil moisture, ice + liquid [m^3 m-3]
!    tah                        - NoahMP36 canopy air temperature [K]
!    eah                        - NoahMP36 canopy air vapor pressure [Pa]
!    fwet                       - NoahMP36 wetted or snowed fraction of canopy [-]
!    canliq                     - NoahMP36 intercepted liquid water [mm]
!    canice                     - NoahMP36 intercepted ice mass [mm]
!    tv                         - NoahMP36 vegetation temperature [K]
!    tg                         - NoahMP36 ground temperature (skin temperature) [K]
!    qsnow                      - NoahMP36 snowfall on the ground [mm s-1]
!    isnow                      - NoahMP36 actual number of snow layers [-]
!    zss                        - NoahMP36 snow/soil layer-bottom depth from snow surface [m]
!    snowh                      - NoahMP36 snow height [m]
!    sneqv                      - NoahMP36 snow water equivalent [mm]
!    snowice                    - NoahMP36 snow-layer ice [mm]
!    snowliq                    - NoahMP36 snow-layer liquid water [mm]
!    zwt                        - NoahMP36 depth to water table [m]
!    wa                         - NoahMP36 water storage in aquifer [mm]
!    wt                         - NoahMP36 water in aquifer and saturated soil [mm]
!    wslake                     - NoahMP36 lake water storage [mm]
!    lfmass                     - NoahMP36 leaf mass [g/m2]
!    rtmass                     - NoahMP36 mass of fine roots [g/m2]
!    stmass                     - NoahMP36 stem mass [g/m2]
!    wood                       - NoahMP36 mass of wood including woody roots [g/m2]
!    stblcp                     - NoahMP36 stable carbon in deep soil [g/m2]
!    fastcp                     - NoahMP36 short-lived carbon in shallow soil [g/m2]
!    lai                        - NoahMP36 leaf area index [-]
!    sai                        - NoahMP36 stem area index [-]
!    cm                         - NoahMP36 momentum drag coefficient [s/m]
!    ch                         - NoahMP36 sensible heat exchange coefficient [s/m]
!    tauss                      - NoahMP36 snow aging term [-]
!    smcwtd                     - NoahMP36 soil water content between bottom of the soil and water table [m^3 m-3]
!    deeprech                   - NoahMP36 recharge to or from the water table when deep [m]
!    rech                       - NoahMP36 recharge to or from the water table when shallow [m]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!      reads a variable from the restart file
!   \item[NoahMP36\_coldstart](\ref{NoahMP36_coldstart}) \newline
!      initializes the NoahMP36 state variables
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
    !WN
    character(len=LIS_CONST_PATH_LEN) :: filen
    integer           :: yr,mo,da,hr,mn,ss,doy
    real*8            :: time
    real              :: gmt
    real              :: ts


    do n=1, LIS_rc%nnest
        wformat = trim(NOAHMP36_struc(n)%rformat)
        ! coldstart
        if(LIS_rc%startcode .eq. "coldstart") then
            call NoahMP36_coldstart(LIS_rc%lsm_index)
        ! restart
        elseif(LIS_rc%startcode .eq. "restart") then
        !WN ---create restart filename based on timewindow for EnKS
                if(LIS_rc%runmode.eq."ensemble smoother") then
                  if(LIS_rc%iterationId(n).gt.1) then
                    if(NOAHMP36_struc(n)%rstInterval.eq.2592000) then
                     !create the restart filename based on the timewindow
                     ! start time
                      call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                           dd=da,calendar=LIS_calendar,rc=status)
                      hr = 0
                      mn = 0
                      ss = 0
                      call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss, &
                           (-1)*LIS_rc%ts)
                    else
                      call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                           dd=da,calendar=LIS_calendar,rc=status)
                      hr = 0
                      mn = 0
                      ss = 0
                    endif

                    call LIS_create_restart_filename(n,filen,'SURFACEMODEL', &
                         'NOAHMP36', &
                         yr,mo,da,hr,mn,ss, wformat=wformat)
                    NOAHMP36_struc(n)%rfile = filen
                  endif
                endif

            allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            ! check the existance of restart file
            inquire(file=NOAHMP36_struc(n)%rfile, exist=file_exists)
            If (.not. file_exists) then
                write(LIS_logunit,*) "NoahMP36 restart file ", &
                     NOAHMP36_struc(n)%rfile," does not exist "
                write(LIS_logunit,*) "Program stopping ..."
                call LIS_endrun
            endif
            write(LIS_logunit,*) "NoahMP36 restart file used: ", &
                 NOAHMP36_struc(n)%rfile

            ! open restart file
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn, file=NOAHMP36_struc(n)%rfile, form="unformatted")
                read(ftn) nc, nr, npatch  !time, veg class, no. tiles

                ! check for grid space conflict
                if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
                    write(LIS_logunit,*) NOAHMP36_struc(n)%rfile, &
                         "grid space mismatch - NoahMP36 halted"
                    call LIS_endrun
                endif

                if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
                    write(LIS_logunit,*) &
                         "restart tile space mismatch, halting..."
                    call LIS_endrun
                endif
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_open(path=NOAHMP36_struc(n)%rfile, &
                                   mode=NF90_NOWRITE, ncid=ftn)
                call LIS_verify(status, &
                     "Error opening file "//NOAHMP36_struc(n)%rfile)
#endif
            endif

            ! read: snow albedo at last time step
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%albold, &
                                     varname="ALBOLD", wformat=wformat)

            ! read: snow mass at the last time step
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%sneqvo, &
                                     varname="SNEQVO", wformat=wformat)

            ! read: snow/soil temperature
            do l=1, NOAHMP36_struc(n)%nsoil + NOAHMP36_struc(n)%nsnow
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="SSTC", &
                     dim=l, &
                     vlevels = NOAHMP36_struc(n)%nsoil + &
                               NOAHMP36_struc(n)%nsnow, &
                     wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP36_struc(n)%noahmp36(t)%sstc(l) = tmptilen(t)
                enddo
            enddo

            ! read: volumetric liquid soil moisture
            do l=1, NOAHMP36_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="SH2O", &
                     dim=l, vlevels = NOAHMP36_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP36_struc(n)%noahmp36(t)%sh2o(l) = tmptilen(t)
                enddo
            enddo

            ! read: volumetric soil moisture, ice + liquid
            do l=1, NOAHMP36_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="SMC", &
                     dim=l, vlevels = NOAHMP36_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                   NOAHMP36_struc(n)%noahmp36(t)%smc(l) = tmptilen(t)
                enddo
            enddo

            ! read: canopy air temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%tah, &
                                     varname="TAH", wformat=wformat)

            ! read: canopy air vapor pressure
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%eah, &
                 varname="EAH", wformat=wformat)

            ! read: wetted or snowed fraction of canopy
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%fwet, &
                 varname="FWET", wformat=wformat)

            ! read: intercepted liquid water
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%canliq, &
                 varname="CANLIQ", wformat=wformat)

            ! read: intercepted ice mass
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%canice, &
                 varname="CANICE", wformat=wformat)

            ! read: vegetation temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%tv, &
                 varname="TV", wformat=wformat)

            ! read: ground temperature (skin temperature)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%tg, &
                 varname="TG", wformat=wformat)

            ! read: snowfall on the ground
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%qsnow, &
                 varname="QSNOW", wformat=wformat)

            ! read: actual number of snow layers
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%isnow, &
                 varname="ISNOW", wformat=wformat)

            ! read: snow/soil layer-bottom depth from snow surface
            do l=1, NOAHMP36_struc(n)%nsoil + NOAHMP36_struc(n)%nsnow
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="ZSS", &
                     dim=l, &
                     vlevels = NOAHMP36_struc(n)%nsoil + &
                               NOAHMP36_struc(n)%nsnow, &
                               wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP36_struc(n)%noahmp36(t)%zss(l) = tmptilen(t)
                enddo
            enddo

            ! read: snow height
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%snowh, &
                 varname="SNOWH", wformat=wformat)

            ! read: snow water equivalent
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%sneqv, &
                 varname="SNEQV", wformat=wformat)

            ! read: snow-layer ice
            do l=1, NOAHMP36_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="SNOWICE", &
                     dim=l, vlevels = NOAHMP36_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP36_struc(n)%noahmp36(t)%snowice(l) = tmptilen(t)
                enddo
            enddo

            ! read: snow-layer liquid water
            do l=1, NOAHMP36_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="SNOWLIQ", &
                     dim=l, vlevels = NOAHMP36_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP36_struc(n)%noahmp36(t)%snowliq(l) = tmptilen(t)
                enddo
            enddo

            ! read: depth to water table
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%zwt, &
                 varname="ZWT", wformat=wformat)

            ! read: water storage in aquifer
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%wa, &
                 varname="WA", wformat=wformat)

            ! read: water in aquifer and saturated soil
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%wt, &
                 varname="WT", wformat=wformat)

            ! read: lake water storage
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%wslake, &
                 varname="WSLAKE", wformat=wformat)

            ! read: leaf mass
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%lfmass, &
                 varname="LFMASS", wformat=wformat)

            ! read: mass of fine roots
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%rtmass, &
                 varname="RTMASS", wformat=wformat)

            ! read: stem mass
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%stmass, &
                 varname="STMASS", wformat=wformat)

            ! read: mass of wood including woody roots
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%wood, &
                 varname="WOOD", wformat=wformat)

            ! read: stable carbon in deep soil
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%stblcp, &
                 varname="STBLCP", wformat=wformat)

            ! read: short-lived carbon in shallow soil
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%fastcp, &
                 varname="FASTCP", wformat=wformat)

            ! read: leaf area index
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%lai, &
                 varname="LAI", wformat=wformat)

            ! read: stem area index
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%sai, &
                 varname="SAI", wformat=wformat)

            ! read: momentum drag coefficient
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%cm, &
                 varname="CM", wformat=wformat)

            ! read: sensible heat exchange coefficient
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%ch, &
                 varname="CH", wformat=wformat)

            ! read: snow aging term
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%tauss, &
                 varname="TAUSS", wformat=wformat)

            ! read: soil water content between bottom of the soil and water
            ! table
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%smcwtd, &
                 varname="SMCWTD", wformat=wformat)

            ! read: recharge to or from the water table when deep
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%deeprech, &
                 varname="DEEPRECH", wformat=wformat)

            ! read: recharge to or from the water table when shallow
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%rech, &
                 varname="RECH", wformat=wformat)

            ! read: reference height for air temperature and humidity
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, &
                 NOAHMP36_struc(n)%noahmp36%zlvl, &
                 varname="ZLVL", wformat=wformat)
            ! close restart file
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, &
                     "Error in nf90_close in NoahMP36_readrst")
#endif
            endif
            deallocate(tmptilen)
        endif
    enddo
end subroutine NoahMP36_readrst

