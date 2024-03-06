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
! !ROUTINE: Ac71_readrst
! \label{Ac71_readrst}
!
! !REVISION HISTORY:
!  06 MAR 2024; Louise Busschaert, initial implementation
!
! !INTERFACE:
subroutine Ac71_readrst()
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify
    use Ac71_lsmMod
    !WN
    use ESMF
    use LIS_fileIOMod
    use LIS_timeMgrMod
    !-----------------
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for Ac71.  This
!  includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the Ac71
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    albold                     - Ac71 snow albedo at last time step [-]
!    sneqvo                     - Ac71 snow mass at the last time step [mm]
!    sstc                       - Ac71 snow/soil temperature [K]
!    sh2o                       - Ac71 volumetric liquid soil moisture [m^3 m-3]
!    smc                        - Ac71 volumetric soil moisture, ice + liquid [m^3 m-3]
!    tah                        - Ac71 canopy air temperature [K]
!    eah                        - Ac71 canopy air vapor pressure [Pa]
!    fwet                       - Ac71 wetted or snowed fraction of canopy [-]
!    canliq                     - Ac71 intercepted liquid water [mm]
!    canice                     - Ac71 intercepted ice mass [mm]
!    tv                         - Ac71 vegetation temperature [K]
!    tg                         - Ac71 ground temperature (skin temperature) [K]
!    qsnow                      - Ac71 snowfall on the ground [mm s-1]
!    isnow                      - Ac71 actual number of snow layers [-]
!    zss                        - Ac71 snow/soil layer-bottom depth from snow surface [m]
!    snowh                      - Ac71 snow height [m]
!    sneqv                      - Ac71 snow water equivalent [mm]
!    snowice                    - Ac71 snow-layer ice [mm]
!    snowliq                    - Ac71 snow-layer liquid water [mm]
!    zwt                        - Ac71 depth to water table [m]
!    wa                         - Ac71 water storage in aquifer [mm]
!    wt                         - Ac71 water in aquifer and saturated soil [mm]
!    wslake                     - Ac71 lake water storage [mm]
!    lfmass                     - Ac71 leaf mass [g/m2]
!    rtmass                     - Ac71 mass of fine roots [g/m2]
!    stmass                     - Ac71 stem mass [g/m2]
!    wood                       - Ac71 mass of wood including woody roots [g/m2]
!    stblcp                     - Ac71 stable carbon in deep soil [g/m2]
!    fastcp                     - Ac71 short-lived carbon in shallow soil [g/m2]
!    lai                        - Ac71 leaf area index [-]
!    sai                        - Ac71 stem area index [-]
!    cm                         - Ac71 momentum drag coefficient [s/m]
!    ch                         - Ac71 sensible heat exchange coefficient [s/m]
!    tauss                      - Ac71 snow aging term [-]
!    smcwtd                     - Ac71 soil water content between bottom of the soil and water table [m^3 m-3]
!    deeprech                   - Ac71 recharge to or from the water table when deep [m]
!    rech                       - Ac71 recharge to or from the water table when shallow [m]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!      reads a variable from the restart file
!   \item[Ac71\_coldstart](\ref{Ac71_coldstart}) \newline
!      initializes the Ac71 state variables
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
    character*100     :: filen
    integer           :: yr,mo,da,hr,mn,ss,doy
    real*8            :: time
    real              :: gmt
    real              :: ts


    do n=1, LIS_rc%nnest
        wformat = trim(AC71_struc(n)%rformat)
        ! coldstart
        if(LIS_rc%startcode .eq. "coldstart") then
            call Ac71_coldstart(LIS_rc%lsm_index)
        ! restart
        elseif(LIS_rc%startcode .eq. "restart") then
        !WN ---create restart filename based on timewindow for EnKS
                if(LIS_rc%runmode.eq."ensemble smoother") then
                  if(LIS_rc%iterationId(n).gt.1) then
                    if(AC71_struc(n)%rstInterval.eq.2592000) then
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
                         'AC71', &
                         yr,mo,da,hr,mn,ss, wformat=wformat)
                    AC71_struc(n)%rfile = filen
                  endif
                endif

            allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            ! check the existance of restart file
            inquire(file=AC71_struc(n)%rfile, exist=file_exists)
            If (.not. file_exists) then
                write(LIS_logunit,*) "Ac71 restart file ", &
                     AC71_struc(n)%rfile," does not exist "
                write(LIS_logunit,*) "Program stopping ..."
                call LIS_endrun
            endif
            write(LIS_logunit,*) "Ac71 restart file used: ", &
                 AC71_struc(n)%rfile

            ! open restart file
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn, file=AC71_struc(n)%rfile, form="unformatted")
                read(ftn) nc, nr, npatch  !time, veg class, no. tiles

                ! check for grid space conflict
                if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
                    write(LIS_logunit,*) AC71_struc(n)%rfile, &
                         "grid space mismatch - Ac71 halted"
                    call LIS_endrun
                endif

                if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
                    write(LIS_logunit,*) &
                         "restart tile space mismatch, halting..."
                    call LIS_endrun
                endif
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_open(path=AC71_struc(n)%rfile, &
                                   mode=NF90_NOWRITE, ncid=ftn)
                call LIS_verify(status, &
                     "Error opening file "//AC71_struc(n)%rfile)
#endif
            endif

            ! read: volumetric soil moisture, ice + liquid
            do l=1, AC71_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="SMC", &
                     dim=l, vlevels = AC71_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                   AC71_struc(n)%ac71(t)%smc(l) = tmptilen(t)
                enddo
            enddo
            ! close restart file
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, &
                     "Error in nf90_close in Ac71_readrst")
#endif
            endif
            deallocate(tmptilen)
        endif
    enddo
end subroutine Ac71_readrst

