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

! AC modules
    use ac_global,  only: IrriMode_Generate, &
                          IrriMode_Inet, &
                          IrriMode_Manual, &
                          IrriMode_NoIrri, &
                          undef_int
    use ac_run,     only: fIrri_open, fIrri_close, &
                          fIrri_read, fIrri_eof
!
! !DESCRIPTION:
!  This program reads restart files for Ac71.  This
!  includes all relevant AC71 variables used to restart.
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
    character*100     :: filen
    integer           :: yr,mo,da,hr,mn,ss,doy
    real*8            :: time
    real              :: gmt
    real              :: ts

    !LB for irrigation rst
    integer           :: i, FromDay_temp, TimeInfo_temp, &
                         DepthInfo_temp, DNr
    character*200     :: TempStr
    real              :: IrriEcw_temp


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

            ! Get irrigation file line number from DayNri
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                if(AC71_struc(n)%ac71(t)%IrriMode.ne.IrriMode_NoIrri) then
                    ! For IrriMode_Generate and IrriMode_Manual
                    if(AC71_struc(n)%ac71(t)%IrriMode.eq.IrriMode_Generate) then
                        ! Initial LN
                        if(AC71_struc(n)%ac71(t)%IrriInfoRecord1%NoMoreInfo)then
                            AC71_struc(n)%ac71(t)%irri_lnr = 11
                        else
                            AC71_struc(n)%ac71(t)%irri_lnr = 12
                            ! re-open irrigation file and read first lines
                            call fIrri_open(trim(AC71_struc(n)%PathNameSimul)&
                                            //trim(AC71_struc(n)%Irrigation_Filename), 'r')
                            do i=1,AC71_struc(n)%ac71(t)%irri_lnr
                                TempStr = fIrri_read()
                            enddo
                            ! Check if we passed the first record
                            do while ((AC71_struc(n)%ac71(t)%daynri-AC71_struc(n)%ac71(t)%Crop%Day1+1)&
                            .gt.AC71_struc(n)%ac71(t)%IrriInfoRecord1%ToDay+1) ! +1 because let the main read
                                ! Read next record
                                TempStr = fIrri_read()
                                ! Extract info (copied and adpated from run.f90, GetIrriParam)
                                if (fIrri_eof()) then
                                    AC71_struc(n)%ac71(t)%IrriInfoRecord1%ToDay = &
                                            AC71_struc(n)%ac71(t)%Crop%DayN - &
                                            AC71_struc(n)%ac71(t)%Crop%Day1 + 1
                                else
                                    AC71_struc(n)%ac71(t)%IrriInfoRecord2%NoMoreInfo = .false.
                                    read(TempStr,*) FromDay_temp, &
                                    TimeInfo_temp, DepthInfo_temp
                                    AC71_struc(n)%ac71(t)%IrriInfoRecord2%FromDay = FromDay_temp
                                    AC71_struc(n)%ac71(t)%IrriInfoRecord2%TimeInfo = TimeInfo_temp
                                    AC71_struc(n)%ac71(t)%IrriInfoRecord2%DepthInfo = DepthInfo_temp
                                    AC71_struc(n)%ac71(t)%Simulation%IrriEcw = IrriEcw_temp
                                end if
                                AC71_struc(n)%ac71(t)%irri_lnr = AC71_struc(n)%ac71(t)%irri_lnr + 1 ! extra record has been read
                            enddo
                        endif
                    elseif(AC71_struc(n)%ac71(t)%IrriMode.eq.IrriMode_Manual) then
                        ! Initial LN
                        if(AC71_struc(n)%ac71(t)%IrriInfoRecord1%NoMoreInfo)then
                            AC71_struc(n)%ac71(t)%irri_lnr = 9
                        else
                            AC71_struc(n)%ac71(t)%irri_lnr = 10
                            ! re-open irrigation file and read the firts lines
                            call fIrri_open(trim(AC71_struc(n)%PathNameSimul)&
                                            //trim(AC71_struc(n)%Irrigation_Filename), 'r')
                            do i=1,AC71_struc(n)%ac71(t)%irri_lnr
                                TempStr = fIrri_read()
                            enddo
                            ! Check if we passed the first record
                            ! Check start date of schedule
                            if(AC71_struc(n)%ac71(t)%IrriFirstDayNr.eq.undef_int)then
                                DNr = AC71_struc(n)%ac71(t)%daynri &
                                    - AC71_struc(n)%ac71(t)%Crop%Day1 + 1
                            else
                                DNr = AC71_struc(n)%ac71(t)%daynri &
                                    - AC71_struc(n)%ac71(t)%IrriFirstDayNr + 1
                            endif
                            do while (DNr.lt.AC71_struc(n)%ac71(t)%IrriInfoRecord1%TimeInfo)
                                ! Read next record
                                TempStr = fIrri_read()
                                ! Extract info (copied and adpated from run.f90, IrriManual)
                                if (fIrri_eof()) then
                                    AC71_struc(n)%ac71(t)%IrriInfoRecord1%NoMoreInfo = .true.
                                else
                                    AC71_struc(n)%ac71(t)%IrriInfoRecord1%NoMoreInfo = .false.
                                    read(TempStr,*) TimeInfo_temp, DepthInfo_temp
                                    AC71_struc(n)%ac71(t)%IrriInfoRecord2%TimeInfo = TimeInfo_temp
                                    AC71_struc(n)%ac71(t)%IrriInfoRecord2%DepthInfo = DepthInfo_temp
                                    AC71_struc(n)%ac71(t)%Simulation%IrriEcw = IrriEcw_temp
                                endif
                                AC71_struc(n)%ac71(t)%irri_lnr = AC71_struc(n)%ac71(t)%irri_lnr + 1 ! extra record has been read
                            enddo
                        endif
                    elseif(AC71_struc(n)%ac71(t)%IrriMode.eq.IrriMode_Inet) then
                        AC71_struc(n)%ac71(t)%irri_lnr = 0 !Inet mode
                    endif
                    ! end irrigation method
                    call fIrri_close() ! close irrigation file
                endif
                ! end check if irrigated

                ! For restart: check if we reached end of simul period
                if (AC71_struc(n)%ac71(t)%DayNri &
                   .eq. AC71_struc(n)%ac71(t)%Simulation%ToDayNr)  then
                    AC71_struc(n)%ac71(t)%irun = 2 ! Means that we need to start a new sim
                    AC71_struc(n)%ac71(t)%InitializeRun = 1
                endif
                
            enddo
            ! end tile
        endif
    enddo
end subroutine Ac71_readrst

