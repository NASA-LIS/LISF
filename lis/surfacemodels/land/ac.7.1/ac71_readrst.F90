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
    integer, allocatable :: tmptilen_int(:)
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
            allocate(tmptilen_int(LIS_rc%npatch(n, LIS_rc%lsm_index)))
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

            ! read: volumetric soil moisture
            do l=1, AC71_struc(n)%max_No_Compartments
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="SMC", &
                     dim=l, vlevels = AC71_struc(n)%max_No_Compartments, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                   AC71_struc(n)%ac71(t)%smc(l) = tmptilen(t)
                   ! Also adjust theta ini
                    !AC71_struc(n)%ac71(t)%Simulation%thetaini(l) = tmptilen(t) !will also be needed for salinity
                enddo
            enddo

            !! From Compartment

            ! read: Compartment_DayAnaero
            do l=1, AC71_struc(n)%max_No_Compartments
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                     varname="Compartment_DayAnaero", &
                     dim=l, vlevels = AC71_struc(n)%max_No_Compartments, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                   AC71_struc(n)%ac71(t)%Compartment(l)%DayAnaero = tmptilen_int(t)
                enddo
            enddo
                                    
            ! read: Compartment_fluxout
            do l=1, AC71_struc(n)%max_No_Compartments
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="Compartment_fluxout", &
                     dim=l, vlevels = AC71_struc(n)%max_No_Compartments, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                   AC71_struc(n)%ac71(t)%Compartment(l)%fluxout = tmptilen(t)
                enddo
            enddo


            ! read: Compartment_Smax
            do l=1, AC71_struc(n)%max_No_Compartments
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="Compartment_Smax", &
                     dim=l, vlevels = AC71_struc(n)%max_No_Compartments, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                   AC71_struc(n)%ac71(t)%Compartment(l)%Smax = tmptilen(t)
                enddo
            enddo

            ! read: Compartment_FCadj
            do l=1, AC71_struc(n)%max_No_Compartments
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="Compartment_FCadj", &
                     dim=l, vlevels = AC71_struc(n)%max_No_Compartments, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                   AC71_struc(n)%ac71(t)%Compartment(l)%FCadj = tmptilen(t)
                enddo
            enddo

            ! read: Compartment_WFactor
            do l=1, AC71_struc(n)%max_No_Compartments
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                     varname="Compartment_WFactor", &
                     dim=l, vlevels = AC71_struc(n)%max_No_Compartments, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                   AC71_struc(n)%ac71(t)%Compartment(l)%WFactor = tmptilen(t)
                enddo
            enddo

            ! read: Compartment_DayAnaero
            do l=1, AC71_struc(n)%max_No_Compartments
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                     varname="Compartment_DayAnaero", &
                     dim=l, vlevels = AC71_struc(n)%max_No_Compartments, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                   AC71_struc(n)%ac71(t)%Compartment(l)%DayAnaero = tmptilen_int(t)
                enddo
            enddo



            !! reals
            ! read: alfaHI
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%alfaHI, &
                                    varname="alfaHI", wformat=wformat)

            ! read: alfaHIAdj
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%alfaHIAdj, &
                                    varname="alfaHIAdj", wformat=wformat)

            ! read: Bin
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%Bin, &
                                    varname="Bin", wformat=wformat)

            ! read: Bout
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%Bout, &
                                    varname="Bout", wformat=wformat)

            ! read: CCiActual
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%CCiActual, &
                                    varname="CCiActual", wformat=wformat)

            ! read: CCiActualWeedInfested
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%CCiActualWeedInfested, &
                                    varname="CCiActualWeedInfested", wformat=wformat)

            ! read: CCiPrev
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%CCiPrev, &
                                    varname="CCiPrev", wformat=wformat)

            ! read: CCiTopEarlySen
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%CCiTopEarlySen, &
                                    varname="CCiTopEarlySen", wformat=wformat)

            ! read: CCxWitheredTpotNoS
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%CCxWitheredTpotNoS, &
                                    varname="CCxWitheredTpotNoS", wformat=wformat)

            ! read: DayFraction
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%DayFraction, &
                                    varname="DayFraction", wformat=wformat)

            ! read: ECstorage
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%ECstorage, &
                                    varname="ECstorage", wformat=wformat)

            ! read: HItimesAT
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%HItimesAT, &
                                    varname="HItimesAT", wformat=wformat)

            ! read: HItimesAT1
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%HItimesAT1, &
                                    varname="HItimesAT1", wformat=wformat)

            ! read: HItimesAT2
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%HItimesAT2, &
                                    varname="HItimesAT2", wformat=wformat)

            ! read: HItimesBEF
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%HItimesBEF, &
                                    varname="HItimesBEF", wformat=wformat)

            ! read: RootZoneWC_Actual
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_Actual, &
                                    varname="RootZoneWC_Actual", wformat=wformat)

            ! read: RootZoneWC_FC
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_FC, &
                                    varname="RootZoneWC_FC", wformat=wformat)

            ! read: RootZoneWC_Leaf
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_Leaf, &
                                    varname="RootZoneWC_Leaf", wformat=wformat)

            ! read: RootZoneWC_SAT
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_SAT, &
                                    varname="RootZoneWC_SAT", wformat=wformat)

            ! read: RootZoneWC_Sen
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_Sen, &
                                    varname="RootZoneWC_Sen", wformat=wformat)

            ! read: RootZoneWC_Thresh
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_Thresh, &
                                    varname="RootZoneWC_Thresh", wformat=wformat)

            ! read: RootZoneWC_WP
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_WP, &
                                    varname="RootZoneWC_WP", wformat=wformat)

            ! read: RootZoneWC_ZtopAct
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_ZtopAct, &
                                    varname="RootZoneWC_ZtopAct", wformat=wformat)

            ! read: RootZoneWC_ZtopFC
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_ZtopFC, &
                                    varname="RootZoneWC_ZtopFC", wformat=wformat)

            ! read: RootZoneWC_ZtopThresh
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_ZtopThresh, &
                                    varname="RootZoneWC_ZtopThresh", wformat=wformat)

            ! read: RootZoneWC_ZtopWP
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%RootZoneWC_ZtopWP, &
                                    varname="RootZoneWC_ZtopWP", wformat=wformat)

            ! read: ScorAT1
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%ScorAT1, &
                                    varname="ScorAT1", wformat=wformat)

            ! read: ScorAT2
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%ScorAT2, &
                                    varname="ScorAT2", wformat=wformat)

            ! read: StressLeaf
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%StressLeaf, &
                                    varname="StressLeaf", wformat=wformat)

            ! read: StressSenescence
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%StressSenescence, &
                                    varname="StressSenescence", wformat=wformat)

            ! read: SumGDDcuts
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%SumGDDcuts, &
                                    varname="SumGDDcuts", wformat=wformat)

            ! read: SumKci
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%SumKci, &
                                    varname="SumKci", wformat=wformat)

            ! read: SumKcTopStress
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%SumKcTopStress, &
                                    varname="SumKcTopStress", wformat=wformat)

            ! read: SurfaceStorage
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%SurfaceStorage, &
                                    varname="SurfaceStorage", wformat=wformat)

            ! read: Tact
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%Tact, &
                                    varname="Tact", wformat=wformat)

            ! read: TactWeedInfested
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%TactWeedInfested, &
                                    varname="TactWeedInfested", wformat=wformat)

            ! read: Tadj
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%Tadj, &
                                    varname="Tadj", wformat=wformat)

            ! read: TimeSenescence
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%TimeSenescence, &
                                    varname="TimeSenescence", wformat=wformat)

            ! read: Tpot
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%Tpot, &
                                    varname="Tpot", wformat=wformat)

            ! read: WeedRCi
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%WeedRCi, &
                                    varname="WeedRCi", wformat=wformat)
            ! read: WPi
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%WPi, &
                                    varname="WPi", wformat=wformat)

            ! read: Ziprev
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%Ziprev, &
                                    varname="Ziprev", wformat=wformat)

            
            !! integers
            ! read: DayNri
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%DayNri, &
                                    varname="DayNri", wformat=wformat)

            ! read: DaySubmerged
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%DaySubmerged, &
                                    varname="DaySubmerged", wformat=wformat)

            ! read: PreviousStressLevel
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%PreviousStressLevel, &
                                    varname="PreviousStressLevel", wformat=wformat)

            ! read: StressSFadjNEW
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%StressSFadjNEW, &
                                    varname="StressSFadjNEW", wformat=wformat)

            ! read: SumInterval
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%SumInterval, &
                                    varname="SumInterval", wformat=wformat)

            !! logical
            !read: NoMoreCrop
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AC71_struc(n)%ac71%NoMoreCrop, &
                                    varname="NoMoreCrop", wformat=wformat)

            !! From StressTot
            ! read: StressTot_Salt
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="StressTot_Salt", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%StressTot%Salt = tmptilen(t)
            enddo

            ! read: StressTot_Temp
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="StressTot_Temp", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%StressTot%Temp = tmptilen(t)
            enddo

            ! read: StressTot_Exp
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="StressTot_Exp", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%StressTot%Exp = tmptilen(t)
            enddo

            ! read: StressTot_Sto
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="StressTot_Sto", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%StressTot%Sto = tmptilen(t)
            enddo

            ! read: StressTot_Weed
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="StressTot_Weed", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%StressTot%Weed = tmptilen(t)
            enddo

            ! read: StressTot_NrD
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                    varname="StressTot_NrD", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%StressTot%NrD = tmptilen_int(t)
            enddo

            !! From SumWaBal
            ! read: SumWaBal_Biomass
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="SumWaBal_Biomass", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%SumWaBal%Biomass = tmptilen(t)
            enddo

            ! read: SumWaBal_BiomassPot
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="SumWaBal_BiomassPot", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%SumWaBal%BiomassPot = tmptilen(t)
            enddo

            ! read: SumWaBal_BiomassTot
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="SumWaBal_BiomassTot", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%SumWaBal%BiomassTot = tmptilen(t)
            enddo

            ! read: SumWaBal_BiomassUnlim
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="SumWaBal_BiomassUnlim", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%SumWaBal%BiomassUnlim = tmptilen(t)
            enddo

            ! read: SumWaBal_YieldPart
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="SumWaBal_YieldPart", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%SumWaBal%YieldPart = tmptilen(t)
            enddo

            !! From Simulation%EffectStress
            ! read: Simulation_EffectStress_CDecline
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="Simulation_EffectStress_CDecline", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Simulation%EffectStress%CDecline = tmptilen(t)
            enddo

            ! read: Simulation_DayAnaero
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                    varname="Simulation_DayAnaero", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Simulation%DayAnaero = tmptilen_int(t)
            enddo

            ! read: Simulation_EffectStress_RedCGC
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                    varname="Simulation_EffectStress_RedCGC", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Simulation%EffectStress%RedCGC = tmptilen_int(t)
            enddo

            ! read: Simulation_EffectStress_RedCCx
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                    varname="Simulation_EffectStress_RedCCx", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Simulation%EffectStress%RedCCx = tmptilen_int(t)
            enddo

            ! read: Simulation_EffectStress_RedWP
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                    varname="Simulation_EffectStress_RedWP", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Simulation%EffectStress%RedWP = tmptilen_int(t)
            enddo

            ! read: Simulation_EffectStress_RedKsSto
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                    varname="Simulation_EffectStress_RedKsSto", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Simulation%EffectStress%RedKsSto = tmptilen_int(t)
            enddo

            ! read: Simulation_SWCtopsoilconsidered
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                    varname="Simulation_SWCtopSoilConsidered", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Simulation%SWCtopSoilConsidered = tmptilen_int(t)
            enddo

            ! read: Simulation_EvapLimitON
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                    varname="Simulation_EvapLimitON", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Simulation%EvapLimitON = tmptilen_int(t)
            enddo

            ! read: Simulation_EvapStartStg2
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                    varname="Simulation_EvapStartStg2", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Simulation%EvapStartStg2 = tmptilen_int(t)
            enddo

            !! From Simulation
            ! read: Simulation_EvapWCSurf
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="Simulation_EvapWCSurf", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Simulation%EvapWCSurf = tmptilen(t)
            enddo

            !! From Management
            ! read: Management_WeedDeltaRC
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_int, &
                    varname="Management_WeedDeltaRC", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Management%WeedDeltaRC = tmptilen_int(t)
            enddo

            !! From Crop
            ! read: Crop_CCxAdjusted
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="Crop_CCxAdjusted", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Crop%CCxAdjusted = tmptilen(t)
            enddo

            ! read: Crop_CCxWithered
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="Crop_CCxWithered", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Crop%CCxWithered = tmptilen(t)
            enddo

            ! read: Crop_pActStom
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="Crop_pActStom", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Crop%pActStom = tmptilen(t)
            enddo

            ! read: Crop_pLeafAct
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="Crop_pLeafAct", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Crop%pLeafAct = tmptilen(t)
            enddo

            ! read: Crop_pSenAct
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                    varname="Crop_pSenAct", wformat=wformat)
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                AC71_struc(n)%ac71(t)%Crop%pSenAct = tmptilen(t)
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

