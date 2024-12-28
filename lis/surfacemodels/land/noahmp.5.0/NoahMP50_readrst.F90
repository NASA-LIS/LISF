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
! !ROUTINE: NoahMP50_readrst
! \label{NoahMP50_readrst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/25/18: Shugong Wang, Zhuo Wang; initial implementation for LIS 7 and Noah-MP-4.0.1
!   01/08/2021 Bailing Li; implemented code for reading GRACE DA restart file
!   May 2023: Cenlin He, modified for refactored NoahMP v5 and later

! !INTERFACE:
subroutine NoahMP50_readrst()
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify                
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use NoahMP50_lsmMod
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
!    sfcrunoff                  - NoahMP accumulated surface runoff [m]
!    udrrunoff                  - NoahMP accumulated sub-surface runoff [m]
!    smc                        - NoahMP volumtric soil moisture [m3/m3]
!    sh2o                       - NoahMP volumtric liquid soil moisture [m3/m3]
!    tslb                       - NoahMP soil temperature [K]
!    sneqv                      - NoahMP snow water equivalent [mm]
!    snowh                      - NoahMP physical snow depth [m]
!    canwat                     - NoahMP total canopy water + ice [mm]
!    acsnom                     - NoahMP accumulated snow melt leaving pack [-]
!    acsnow                     - NoahMP accumulated snow on grid [mm]
!    isnow                      - NoahMP actual no. of snow layers [-]
!    tv                         - NoahMP vegetation leaf temperature [K]
!    tg                         - NoahMP bulk ground surface temperature [K]
!    canice                     - NoahMP canopy-intercepted ice [mm]
!    canliq                     - NoahMP canopy-intercepted liquid water [mm]
!    eah                        - NoahMP canopy air vapor pressure [Pa]
!    tah                        - NoahMP canopy air temperature [K]
!    cm                         - NoahMP bulk momentum drag coefficient [-]
!    ch                         - NoahMP bulk sensible heat exchange coefficient [-]
!    fwet                       - NoahMP wetted or snowed fraction of canopy [-]
!    sneqvo                     - NoahMP snow mass at last time step [mm h2o]
!    albold                     - NoahMP snow albedo at last time step [-]
!    qsnow                      - NoahMP snowfall on the ground [mm/s]
!    wslake                     - NoahMP lake water storage [mm]
!    zwt                        - NoahMP water table depth [m]
!    wa                         - NoahMP water in the "aquifer" [mm]
!    wt                         - NoahMP water in aquifer and saturated soil [mm]
!    tsno                       - NoahMP snow layer temperature [K]
!    zss                        - NoahMP snow/soil layer depth from snow surface [m]
!    snowice                    - NoahMP snow layer ice [mm]
!    snowliq                    - NoahMP snow layer liquid water [mm]
!    lfmass                     - NoahMP leaf mass [g/m2]
!    rtmass                     - NoahMP mass of fine roots [g/m2]
!    stmass                     - NoahMP stem mass [g/m2]
!    wood                       - NoahMP mass of wood (including woody roots) [g/m2]
!    stblcp                     - NoahMP stable carbon in deep soil [g/m2]
!    fastcp                     - NoahMP short-lived carbon in shallow soil [g/m2]
!    lai                        - NoahMP leaf area index [-]
!    sai                        - NoahMP stem area index [-]
!    tauss                      - NoahMP snow age factor [-]
!    smoiseq                    - NoahMP equilibrium volumetric soil moisture content [m3/m3]
!    smcwtd                     - NoahMP soil moisture content in the layer to the water table when deep [-]
!    deeprech                   - NoahMP recharge to the water table when deep [-]
!    rech                       - NoahMP recharge to the water table (diagnostic) [-]
!    grain                      - NoahMP mass of grain XING [g/m2]
!    gdd                        - NoahMP growing degree days XING (based on 10C) [-]
!    pgs                        - NoahMP growing degree days XING [-]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart})\\
!      reads a variable from the restart file
!   \item[NoahMP50\_coldstart](\ref{NoahMP50_coldstart})\\
!      initializes the NoahMP state variables
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
        wformat = trim(NoahMP50_struc(n)%rformat)
        ! coldstart
        if(LIS_rc%startcode .eq. "coldstart") then  
            call NoahMP50_coldstart(LIS_rc%lsm_index)
        ! restart
        elseif(LIS_rc%startcode .eq. "restart") then
        !---create restart filename based on timewindow for EnKS
                if(LIS_rc%runmode.eq."ensemble smoother") then
                  if(LIS_rc%iterationId(n).gt.1) then
                    if(NoahMP50_struc(n)%rstInterval.eq.2592000) then
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

                    call LIS_create_restart_filename(n,filen,'SURFACEMODEL','NoahMP50', &
                         yr,mo,da,hr,mn,ss, wformat=wformat)
                    NoahMP50_struc(n)%rfile = filen
                  endif
                endif


            allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            ! check the existance of restart file
            inquire(file=NoahMP50_struc(n)%rfile, exist=file_exists)
            If (.not. file_exists) then 
                write(LIS_logunit,*) "[ERR] Noah-MP.5.0 restart file: ", &
                                      trim(NoahMP50_struc(n)%rfile)
                write(LIS_logunit,*) "[ERR] does not exist."
                write(LIS_logunit,*) "[ERR] Program stopping ..."
                call LIS_endrun
            endif
            write(LIS_logunit,*) &
                 "[INFO] Noah-MP.5.0 restart file used: ",trim(NoahMP50_struc(n)%rfile)
        
            ! open restart file
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn, file=NoahMP50_struc(n)%rfile, &
                     form="unformatted")
                read(ftn) nc, nr, npatch  !time, veg class, no. tiles
 
                ! check for grid space conflict
                if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
                   write(LIS_logunit,*) "[ERR]",trim(NoahMP50_struc(n)%rfile)
                   write(LIS_logunit,*) "[ERR] grid space mismatch"
                   write(LIS_logunit,*) "[ERR] Program stopping..."
                   call LIS_endrun
                endif
            
                if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
                   write(LIS_logunit,*) "[ERR]",trim(NoahMP50_struc(n)%rfile)
                   write(LIS_logunit,*) "[ERR] tile space mismatch"
                   write(LIS_logunit,*) "[ERR] Program stopping..."
                   call LIS_endrun
                endif
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_open(path=NoahMP50_struc(n)%rfile, &
                                   mode=NF90_NOWRITE, ncid=ftn)
                call LIS_verify(status, "Error opening file "//NoahMP50_struc(n)%rfile)
#endif
            endif
 
            ! read: accumulated surface runoff
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%sfcrunoff, &
                                     varname="SFCRUNOFF", wformat=wformat)
 
            ! read: accumulated sub-surface runoff
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%udrrunoff, &
                                     varname="UDRRUNOFF", wformat=wformat)
 
            ! read: volumtric soil moisture
            do l=1, NoahMP50_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SMC", &
                                         dim=l, vlevels = NoahMP50_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NoahMP50_struc(n)%noahmp50(t)%smc(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: volumtric liquid soil moisture
            do l=1, NoahMP50_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SH2O", &
                                         dim=l, vlevels = NoahMP50_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NoahMP50_struc(n)%noahmp50(t)%sh2o(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: soil temperature
            do l=1, NoahMP50_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="TSLB", &
                                         dim=l, vlevels = NoahMP50_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NoahMP50_struc(n)%noahmp50(t)%tslb(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: snow water equivalent
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%sneqv, &
                                     varname="SNEQV", wformat=wformat)
 
            ! read: physical snow depth
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%snowh, &
                                     varname="SNOWH", wformat=wformat)
 
            ! read: total canopy water + ice
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%canwat, &
                                     varname="CANWAT", wformat=wformat)
 
            ! read: accumulated snow melt leaving pack
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%acsnom, &
                                     varname="ACSNOM", wformat=wformat)
 
            ! read: accumulated snow on grid
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%acsnow, &
                                     varname="ACSNOW", wformat=wformat)
 
            ! read: actual no. of snow layers
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%isnow, &
                                     varname="ISNOW", wformat=wformat)
 
            ! read: vegetation leaf temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%tv, &
                                     varname="TV", wformat=wformat)
 
            ! read: bulk ground surface temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%tg, &
                                     varname="TG", wformat=wformat)
 
            ! read: canopy-intercepted ice
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%canice, &
                                     varname="CANICE", wformat=wformat)
 
            ! read: canopy-intercepted liquid water
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%canliq, &
                                     varname="CANLIQ", wformat=wformat)
 
            ! read: canopy air vapor pressure
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%eah, &
                                     varname="EAH", wformat=wformat)
 
            ! read: canopy air temperature
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%tah, &
                                     varname="TAH", wformat=wformat)
 
            ! read: bulk momentum drag coefficient
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%cm, &
                                     varname="CM", wformat=wformat)
 
            ! read: bulk sensible heat exchange coefficient
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%ch, &
                                     varname="CH", wformat=wformat)
 
            ! read: wetted or snowed fraction of canopy
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%fwet, &
                                     varname="FWET", wformat=wformat)
 
            ! read: snow mass at last time step
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%sneqvo, &
                                     varname="SNEQVO", wformat=wformat)
 
            ! read: snow albedo at last time step
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%albold, &
                                     varname="ALBOLD", wformat=wformat)
 
            ! read: snowfall on the ground
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%qsnow, &
                                     varname="QSNOW", wformat=wformat)
 
            ! read: lake water storage
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%wslake, &
                                     varname="WSLAKE", wformat=wformat)
 
            ! read: water table depth
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%zwt, &
                                     varname="ZWT", wformat=wformat)
 
            ! read: water in the "aquifer"
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%wa, &
                                     varname="WA", wformat=wformat)
 
            ! read: water in aquifer and saturated soil
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%wt, &
                                     varname="WT", wformat=wformat)
 
            ! read: snow layer temperature
            do l=1, NoahMP50_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="TSNO", &
                                         dim=l, vlevels = NoahMP50_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NoahMP50_struc(n)%noahmp50(t)%tsno(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: snow/soil layer depth from snow surface
            do l=1, NoahMP50_struc(n)%nsnow + NoahMP50_struc(n)%nsoil
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="ZSS", &
                                         dim=l, vlevels = NoahMP50_struc(n)%nsnow + NoahMP50_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NoahMP50_struc(n)%noahmp50(t)%zss(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: snow layer ice
            do l=1, NoahMP50_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWICE", &
                                         dim=l, vlevels = NoahMP50_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NoahMP50_struc(n)%noahmp50(t)%snowice(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: snow layer liquid water
            do l=1, NoahMP50_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWLIQ", &
                                         dim=l, vlevels = NoahMP50_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NoahMP50_struc(n)%noahmp50(t)%snowliq(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: leaf mass
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%lfmass, &
                                     varname="LFMASS", wformat=wformat)
 
            ! read: mass of fine roots
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%rtmass, &
                                     varname="RTMASS", wformat=wformat)
 
            ! read: stem mass
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%stmass, &
                                     varname="STMASS", wformat=wformat)
 
            ! read: mass of wood (including woody roots)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%wood, &
                                     varname="WOOD", wformat=wformat)
 
            ! read: stable carbon in deep soil
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%stblcp, &
                                     varname="STBLCP", wformat=wformat)
 
            ! read: short-lived carbon in shallow soil
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%fastcp, &
                                     varname="FASTCP", wformat=wformat)
 
            ! read: leaf area index
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%lai, &
                                     varname="LAI", wformat=wformat)
 
            ! read: stem area index
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%sai, &
                                     varname="SAI", wformat=wformat)
 
            ! read: snow age factor
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%tauss, &
                                     varname="TAUSS", wformat=wformat)

            if (NoahMP50_struc(n)%runsub_opt == 5) then 
               ! read: equilibrium volumetric soil moisture content
               do l=1, NoahMP50_struc(n)%nsoil ! TODO: check loop
                  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SMOISEQ", &
                                           dim=l, vlevels = NoahMP50_struc(n)%nsoil, wformat=wformat)
                   do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                      NoahMP50_struc(n)%noahmp50(t)%smoiseq(l) = tmptilen(t)
                   enddo
               enddo
 
               ! read: soil moisture content in the layer to the water table when deep
               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%smcwtd, &
                                        varname="SMCWTD", wformat=wformat)
 
               ! read: recharge to the water table when deep
               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%deeprech, &
                                        varname="DEEPRECH", wformat=wformat)
 
               ! read: recharge to the water table (diagnostic)
               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%rech, &
                                        varname="RECH", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%pexp, &
                                        varname="PEXP", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%area, &
                                        varname="AREA", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%qrf, &
                                        varname="QRF", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%qspring, &
                                        varname="QSPRING", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%qslat, &
                                        varname="QSLAT", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%qrfs, &
                                        varname="QRFS", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%qsprings, &
                                        varname="QSPRINGS", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%fdepth, &
                                        varname="FDEPTH", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%rivercond, &
                                        varname="RIVERCOND", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%riverbed, &
                                        varname="RIVERBED", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%eqzwt, &
                                        varname="EQZWT", wformat=wformat)
            endif ! MMF groundwater 

            ! for irrigation
            if (NoahMP50_struc(n)%irr_opt >0) then
               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%irnumsi, &
                                        varname="IRNUMSI", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%irnummi, &
                                        varname="IRNUMMI", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%irnumfi, &
                                        varname="IRNUMFI", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%irwatsi, &
                                        varname="IRWATSI", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%irwatmi, &
                                        varname="IRWATMI", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%irwatfi, &
                                        varname="IRWATFI", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%irsivol, &
                                        varname="IRSIVOL", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%irmivol, &
                                        varname="IRMIVOL", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%irfivol, &
                                        varname="IRFIVOL", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%ireloss, &
                                        varname="IRELOSS", wformat=wformat)

               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%irrsplh, &
                                        varname="IRRSPLH", wformat=wformat)
            endif ! irrigation 

            ! for tile drainage
            if (NoahMP50_struc(n)%tdrn_opt >0) then
               call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%qtdrain, &
                                        varname="QTDRAIN", wformat=wformat)
            endif

            ! for crop 
            ! read: mass of grain XING
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%grain, &
                                     varname="GRAIN", wformat=wformat)
 
            ! read: growing degree days XING (based on 10C)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%gdd, &
                                     varname="GDD", wformat=wformat)
 
            ! read: growing degree days XING
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%pgs, &
                                     varname="PGS", wformat=wformat)

            ! for additional restart variables
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%accssoil, &
                                     varname="ACC_SSOIL", wformat=wformat)

            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%accqinsur, &
                                     varname="ACC_QINSUR", wformat=wformat)

            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%accqseva, &
                                     varname="ACC_QSEVA", wformat=wformat)

            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%accdwater, &
                                     varname="ACC_DWATER", wformat=wformat)

            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%accprcp, &
                                     varname="ACC_PRCP", wformat=wformat)

            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%accecan, &
                                     varname="ACC_ECAN", wformat=wformat)

            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%accetran, &
                                     varname="ACC_ETRAN", wformat=wformat)

            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, NoahMP50_struc(n)%noahmp50%accedir, &
                                     varname="ACC_EDIR", wformat=wformat)

            do l=1, NoahMP50_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="ACC_ETRANI", &
                                         dim=l, vlevels = NoahMP50_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    NoahMP50_struc(n)%noahmp50(t)%accetrani(l) = tmptilen(t)
                enddo
            enddo

            ! close restart file
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, &
                     "Error in nf90_close in NoahMP50_readrst")
#endif
            endif
            deallocate(tmptilen)
        endif    
    enddo
end subroutine NoahMP50_readrst 
