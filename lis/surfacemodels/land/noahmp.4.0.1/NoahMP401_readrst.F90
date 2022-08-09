!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LIS_misc.h"
!BOP
!
! !ROUTINE: NoahMP401_readrst
! \label{NoahMP401_readrst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/25/18: Shugong Wang, Zhuo Wang; initial implementation for LIS 7 and Noah-MP-4.0.1
!   01/08/2021 Bailing Li; implemented code for reading GRACE DA restart file
!
! !INTERFACE:
subroutine NoahMP401_readrst()
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify                
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use NoahMP401_lsmMod
    use ESMF
    use LIS_fileIOMod
    use LIS_timeMgrMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for Noah-MP-4.0.1 LSM.
!  This includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the Noah-MP-4.0.1 
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    sfcrunoff                  - NoahMP401 accumulated surface runoff [m]
!    udrrunoff                  - NoahMP401 accumulated sub-surface runoff [m]
!    smc                        - NoahMP401 volumtric soil moisture [m3/m3]
!    sh2o                       - NoahMP401 volumtric liquid soil moisture [m3/m3]
!    tslb                       - NoahMP401 soil temperature [K]
!    sneqv                      - NoahMP401 snow water equivalent [mm]
!    snowh                      - NoahMP401 physical snow depth [m]
!    canwat                     - NoahMP401 total canopy water + ice [mm]
!    acsnom                     - NoahMP401 accumulated snow melt leaving pack [-]
!    acsnow                     - NoahMP401 accumulated snow on grid [mm]
!    isnow                      - NoahMP401 actual no. of snow layers [-]
!    tv                         - NoahMP401 vegetation leaf temperature [K]
!    tg                         - NoahMP401 bulk ground surface temperature [K]
!    canice                     - NoahMP401 canopy-intercepted ice [mm]
!    canliq                     - NoahMP401 canopy-intercepted liquid water [mm]
!    eah                        - NoahMP401 canopy air vapor pressure [Pa]
!    tah                        - NoahMP401 canopy air temperature [K]
!    cm                         - NoahMP401 bulk momentum drag coefficient [-]
!    ch                         - NoahMP401 bulk sensible heat exchange coefficient [-]
!    fwet                       - NoahMP401 wetted or snowed fraction of canopy [-]
!    sneqvo                     - NoahMP401 snow mass at last time step [mm h2o]
!    albold                     - NoahMP401 snow albedo at last time step [-]
!    qsnow                      - NoahMP401 snowfall on the ground [mm/s]
!    wslake                     - NoahMP401 lake water storage [mm]
!    zwt                        - NoahMP401 water table depth [m]
!    wa                         - NoahMP401 water in the "aquifer" [mm]
!    wt                         - NoahMP401 water in aquifer and saturated soil [mm]
!    tsno                       - NoahMP401 snow layer temperature [K]
!    zss                        - NoahMP401 snow/soil layer depth from snow surface [m]
!    snowice                    - NoahMP401 snow layer ice [mm]
!    snowliq                    - NoahMP401 snow layer liquid water [mm]
!    lfmass                     - NoahMP401 leaf mass [g/m2]
!    rtmass                     - NoahMP401 mass of fine roots [g/m2]
!    stmass                     - NoahMP401 stem mass [g/m2]
!    wood                       - NoahMP401 mass of wood (including woody roots) [g/m2]
!    stblcp                     - NoahMP401 stable carbon in deep soil [g/m2]
!    fastcp                     - NoahMP401 short-lived carbon in shallow soil [g/m2]
!    lai                        - NoahMP401 leaf area index [-]
!    sai                        - NoahMP401 stem area index [-]
!    tauss                      - NoahMP401 snow age factor [-]
!    smoiseq                    - NoahMP401 equilibrium volumetric soil moisture content [m3/m3]
!    smcwtd                     - NoahMP401 soil moisture content in the layer to the water table when deep [-]
!    deeprech                   - NoahMP401 recharge to the water table when deep [-]
!    rech                       - NoahMP401 recharge to the water table (diagnostic) [-]
!    grain                      - NoahMP401 mass of grain XING [g/m2]
!    gdd                        - NoahMP401 growing degree days XING (based on 10C) [-]
!    pgs                        - NoahMP401 growing degree days XING [-]
!    gecros_state               - NoahMP401 optional gecros crop [-]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart})\\
!      reads a variable from the restart file
!   \item[NoahMP401\_coldstart](\ref{NoahMP401_coldstart})\\
!      initializes the NoahMP401 state variables
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
    character(len=LIS_CONST_PATH_LEN) :: filen
    integer           :: yr,mo,da,hr,mn,ss,doy
    real*8            :: time
    real              :: gmt
    real              :: ts

 
    do n=1, LIS_rc%nnest
        wformat = trim(NOAHMP401_struc(n)%rformat)
        ! coldstart
        if(LIS_rc%startcode .eq. "coldstart") then  
            call NoahMP401_coldstart(LIS_rc%lsm_index)
        ! restart
        elseif(LIS_rc%startcode .eq. "restart") then
        !---create restart filename based on timewindow for EnKS
                if(LIS_rc%runmode.eq."ensemble smoother") then
                  if(LIS_rc%iterationId(n).gt.1) then
                    if(NOAHMP401_struc(n)%rstInterval.eq.2592000) then
                     !create the restart filename based on the timewindow start time
                      call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                           dd=da,calendar=LIS_calendar,rc=status)
                      hr = 0
                      mn = 0
                      ss = 0
                      call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,(-1)*LIS_rc%ts)
                    else
                      call ESMF_TimeGet(LIS_twStartTime,yy=yr,mm=mo,&
                           dd=da,calendar=LIS_calendar,rc=status)
                      hr = 0
                      mn = 0
                      ss = 0
                    endif

                    call LIS_create_restart_filename(n,filen,'SURFACEMODEL','NOAHMP401', &
                         yr,mo,da,hr,mn,ss, wformat=wformat)
                    NOAHMP401_struc(n)%rfile = filen
                  endif
                endif


            allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            ! check the existance of restart file
            inquire(file=NOAHMP401_struc(n)%rfile, exist=file_exists)
            If (.not. file_exists) then 
                write(LIS_logunit,*) "[ERR] NoahMP401 restart file: ", &
                                      trim(NOAHMP401_struc(n)%rfile)
                write(LIS_logunit,*) "[ERR] does not exist."
                write(LIS_logunit,*) "[ERR] Program stopping ..."
                call LIS_endrun
            endif
            write(LIS_logunit,*) &
                 "[INFO] Noah-MP.4.0.1 restart file used: ",trim(NOAHMP401_struc(n)%rfile)
        
            ! open restart file
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn, file=NOAHMP401_struc(n)%rfile, &
                     form="unformatted")
                read(ftn) nc, nr, npatch  !time, veg class, no. tiles
 
                ! check for grid space conflict
                if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
                   write(LIS_logunit,*) "[ERR]",trim(NOAHMP401_struc(n)%rfile)
                   write(LIS_logunit,*) "[ERR] grid space mismatch"
                   write(LIS_logunit,*) "[ERR] Program stopping..."
                   call LIS_endrun
                endif
            
                if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
                   write(LIS_logunit,*) "[ERR]",trim(NOAHMP401_struc(n)%rfile)
                   write(LIS_logunit,*) "[ERR] tile space mismatch"
                   write(LIS_logunit,*) "[ERR] Program stopping..."
                   call LIS_endrun
                endif
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_open(path=NOAHMP401_struc(n)%rfile, &
                                   mode=NF90_NOWRITE, ncid=ftn)
                call LIS_verify(status, "Error opening file "//NOAHMP401_struc(n)%rfile)
#endif
            endif
 
            ! read: accumulated surface runoff
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%sfcrunoff, &
                                     varname="SFCRUNOFF", wformat=wformat)
 
            ! read: accumulated sub-surface runoff
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%udrrunoff, &
                                     varname="UDRRUNOFF", wformat=wformat)
 
            ! read: volumtric soil moisture
            do l=1, NOAHMP401_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SMC", &
                                         dim=l, vlevels = NOAHMP401_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP401_struc(n)%noahmp401(t)%smc(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: volumtric liquid soil moisture
            do l=1, NOAHMP401_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SH2O", &
                                         dim=l, vlevels = NOAHMP401_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP401_struc(n)%noahmp401(t)%sh2o(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: soil temperature
            do l=1, NOAHMP401_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="TSLB", &
                                         dim=l, vlevels = NOAHMP401_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP401_struc(n)%noahmp401(t)%tslb(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: snow water equivalent
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%sneqv, &
                                     varname="SNEQV", wformat=wformat)
 
            ! read: physical snow depth
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%snowh, &
                                     varname="SNOWH", wformat=wformat)
 
            ! read: total canopy water + ice
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%canwat, &
                                     varname="CANWAT", wformat=wformat)
 
            ! read: accumulated snow melt leaving pack
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%acsnom, &
                                     varname="ACSNOM", wformat=wformat)
 
            ! read: accumulated snow on grid
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%acsnow, &
                                     varname="ACSNOW", wformat=wformat)
 
            ! read: actual no. of snow layers
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%isnow, &
                                     varname="ISNOW", wformat=wformat)
 
            ! read: vegetation leaf temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%tv, &
                                     varname="TV", wformat=wformat)
 
            ! read: bulk ground surface temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%tg, &
                                     varname="TG", wformat=wformat)
 
            ! read: canopy-intercepted ice
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%canice, &
                                     varname="CANICE", wformat=wformat)
 
            ! read: canopy-intercepted liquid water
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%canliq, &
                                     varname="CANLIQ", wformat=wformat)
 
            ! read: canopy air vapor pressure
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%eah, &
                                     varname="EAH", wformat=wformat)
 
            ! read: canopy air temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%tah, &
                                     varname="TAH", wformat=wformat)
 
            ! read: bulk momentum drag coefficient
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%cm, &
                                     varname="CM", wformat=wformat)
 
            ! read: bulk sensible heat exchange coefficient
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%ch, &
                                     varname="CH", wformat=wformat)
 
            ! read: wetted or snowed fraction of canopy
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%fwet, &
                                     varname="FWET", wformat=wformat)
 
            ! read: snow mass at last time step
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%sneqvo, &
                                     varname="SNEQVO", wformat=wformat)
 
            ! read: snow albedo at last time step
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%albold, &
                                     varname="ALBOLD", wformat=wformat)
 
            ! read: snowfall on the ground
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%qsnow, &
                                     varname="QSNOW", wformat=wformat)
 
            ! read: lake water storage
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%wslake, &
                                     varname="WSLAKE", wformat=wformat)
 
            ! read: water table depth
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%zwt, &
                                     varname="ZWT", wformat=wformat)
 
            ! read: water in the "aquifer"
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%wa, &
                                     varname="WA", wformat=wformat)
 
            ! read: water in aquifer and saturated soil
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%wt, &
                                     varname="WT", wformat=wformat)
 
            ! read: snow layer temperature
            do l=1, NOAHMP401_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="TSNO", &
                                         dim=l, vlevels = NOAHMP401_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP401_struc(n)%noahmp401(t)%tsno(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: snow/soil layer depth from snow surface
            do l=1, NOAHMP401_struc(n)%nsnow + NOAHMP401_struc(n)%nsoil
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="ZSS", &
                                         dim=l, vlevels = NOAHMP401_struc(n)%nsnow + NOAHMP401_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP401_struc(n)%noahmp401(t)%zss(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: snow layer ice
            do l=1, NOAHMP401_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWICE", &
                                         dim=l, vlevels = NOAHMP401_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP401_struc(n)%noahmp401(t)%snowice(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: snow layer liquid water
            do l=1, NOAHMP401_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWLIQ", &
                                         dim=l, vlevels = NOAHMP401_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP401_struc(n)%noahmp401(t)%snowliq(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: leaf mass
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%lfmass, &
                                     varname="LFMASS", wformat=wformat)
 
            ! read: mass of fine roots
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%rtmass, &
                                     varname="RTMASS", wformat=wformat)
 
            ! read: stem mass
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%stmass, &
                                     varname="STMASS", wformat=wformat)
 
            ! read: mass of wood (including woody roots)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%wood, &
                                     varname="WOOD", wformat=wformat)
 
            ! read: stable carbon in deep soil
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%stblcp, &
                                     varname="STBLCP", wformat=wformat)
 
            ! read: short-lived carbon in shallow soil
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%fastcp, &
                                     varname="FASTCP", wformat=wformat)
 
            ! read: leaf area index
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%lai, &
                                     varname="LAI", wformat=wformat)
 
            ! read: stem area index
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%sai, &
                                     varname="SAI", wformat=wformat)
 
            ! read: snow age factor
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%tauss, &
                                     varname="TAUSS", wformat=wformat)
 
            ! read: equilibrium volumetric soil moisture content
            do l=1, NOAHMP401_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SMOISEQ", &
                                         dim=l, vlevels = NOAHMP401_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NOAHMP401_struc(n)%noahmp401(t)%smoiseq(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: soil moisture content in the layer to the water table when deep
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%smcwtd, &
                                     varname="SMCWTD", wformat=wformat)
 
            ! read: recharge to the water table when deep
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%deeprech, &
                                     varname="DEEPRECH", wformat=wformat)
 
            ! read: recharge to the water table (diagnostic)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%rech, &
                                     varname="RECH", wformat=wformat)
 
            ! read: mass of grain XING
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%grain, &
                                     varname="GRAIN", wformat=wformat)
 
            ! read: growing degree days XING (based on 10C)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%gdd, &
                                     varname="GDD", wformat=wformat)
 
            ! read: growing degree days XING
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP401_struc(n)%noahmp401%pgs, &
                                     varname="PGS", wformat=wformat)
 
            ! read: optional gecros crop
!            do l=1, 60 ! TODO: check loop
!                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="GECROS_STATE", &
!                                         dim=l, vlevels = 60, wformat=wformat)
!                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
!                    NOAHMP401_struc(n)%noahmp401(t)%gecros_state(l) = tmptilen(t)
!                enddo
!            enddo
        
            ! close restart file
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, &
                     "Error in nf90_close in NoahMP401_readrst")
#endif
            endif
            deallocate(tmptilen)
        endif    
    enddo
end subroutine NoahMP401_readrst
        
