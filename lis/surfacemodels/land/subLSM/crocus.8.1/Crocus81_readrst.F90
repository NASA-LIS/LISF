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
! !ROUTINE: Crocus81_readrst
! \label{Crocus81_readrst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/18/19: Mahdi Navari, Shugong Wang; initial implementation for LIS 7 and Crocus81
!
! !INTERFACE:
subroutine Crocus81_readrst()
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify                
    use Crocus81_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for Crocus81.  This
!  includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the Crocus81
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    SNOWSWE                    - Crocus81 Snow layer(s) liquid Water Equivalent (SWE:kg m-2) [kg/m2]
!    SNOWRHO                    - Crocus81 Snow layer(s) averaged density (kg/m3) [kg/m3]
!    SNOWHEAT                   - Crocus81 Snow layer(s) heat content (J/m2) [J/m2]
!    SNOWALB                    - Crocus81 snow surface albedo [-]
!    SNOWGRAN1                  - Crocus81 Snow layers grain feature 1 [-]
!    SNOWGRAN2                  - Crocus81 Snow layer grain feature 2 [-]
!    SNOWHIST                   - Crocus81 Snow layer grain historical parameter (only for non dendritic snow) (-) in {0-5} [-]
!    SNOWAGE                    - Crocus81 Age since snowfall (day) [day]
!    SNOWLIQ                    - Crocus81 Snow layer(s) liquid water content (m) [m]
!    SNOWTEMP                   - Crocus81 Snow layer(s) temperature (K) [K]
!    SNOWDZ                     - Crocus81 Snow layer(s) thickness (m) [m]
!    GRNDFLUX                   - Crocus81 Soil/snow interface heat flux (W/m2) [W/m2]
!    SNDRIFT                    - Crocus81 Blowing snow sublimation (kg/m2/s) NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592) [kg/m2/s]
!    RI_n                       - Crocus81 Richardson number (-)  NOTE: RI has not been initialized in CALL_MODEL (If not OMED initalized to undefined in the snow3L_isba.F90) [-]
!    CDSNOW                     - Crocus81 Drag coefficient for momentum over snow (-) [-]
!    USTARSNOW                  - Crocus81 Friction velocity over snow (m/s); [m/s]
!    CHSNOW                     - Crocus81 Drag coefficient for heat over snow  (-) [-]
!    SNOWMAK_dz                 - Crocus81 Snowmaking thickness (m) [m]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart})\\
!      reads a variable from the restart file
!   \item[Crocus81\_coldstart](\ref{Crocus81_coldstart})\\
!      initializes the Crocus81 state variables
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
        wformat = trim(CROCUS81_struc(n)%rformat)
        ! coldstart
        if(LIS_rc%startcode .eq. "coldstart") then  
            call Crocus81_coldstart(LIS_rc%lsm_index)
        ! restart
        elseif(LIS_rc%startcode .eq. "restart") then
            allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            ! check the existance of restart file
            inquire(file=CROCUS81_struc(n)%rfile, exist=file_exists)
            If (.not. file_exists) then 
                write(LIS_logunit,*) "Crocus81 restart file ", CROCUS81_struc(n)%rfile," does not exist "
                write(LIS_logunit,*) "Program stopping ..."
                call LIS_endrun
            endif
            write(LIS_logunit,*) "Crocus81 restart file used: ", CROCUS81_struc(n)%rfile
        
            ! open restart file
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn, file=CROCUS81_struc(n)%rfile, form="unformatted")
                read(ftn) nc, nr, npatch  !time, veg class, no. tiles
 
                ! check for grid space conflict
                if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
                    write(LIS_logunit,*) CROCUS81_struc(n)%rfile, "grid space mismatch - Crocus81 halted"
                    call LIS_endrun
                endif
            
                if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
                    write(LIS_logunit,*) "restart tile space mismatch, halting..."
                    call LIS_endrun
                endif
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_open(path=CROCUS81_struc(n)%rfile, &
                                   mode=NF90_NOWRITE, ncid=ftn)
                call LIS_verify(status, "Error opening file "//CROCUS81_struc(n)%rfile)
#endif
            endif
 
            ! read: Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
            do l=1, CROCUS81_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWSWE", &
                                         dim=l, vlevels = CROCUS81_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    CROCUS81_struc(n)%crocus81(t)%SNOWSWE(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer(s) averaged density (kg/m3)
            do l=1, CROCUS81_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWRHO", &
                                         dim=l, vlevels = CROCUS81_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    CROCUS81_struc(n)%crocus81(t)%SNOWRHO(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer(s) heat content (J/m2)
            do l=1, CROCUS81_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWHEAT", &
                                         dim=l, vlevels = CROCUS81_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    CROCUS81_struc(n)%crocus81(t)%SNOWHEAT(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: snow surface albedo
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%SNOWALB, &
                                     varname="SNOWALB", wformat=wformat)
 
            ! read: Snow layers grain feature 1
            do l=1, CROCUS81_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWGRAN1", &
                                         dim=l, vlevels = CROCUS81_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    CROCUS81_struc(n)%crocus81(t)%SNOWGRAN1(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer grain feature 2
            do l=1, CROCUS81_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWGRAN2", &
                                         dim=l, vlevels = CROCUS81_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer grain historical parameter (only for non dendritic snow) (-) in {0-5}
            do l=1, CROCUS81_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWHIST", &
                                         dim=l, vlevels = CROCUS81_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    CROCUS81_struc(n)%crocus81(t)%SNOWHIST(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Age since snowfall (day)
            do l=1, CROCUS81_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWAGE", &
                                         dim=l, vlevels = CROCUS81_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    CROCUS81_struc(n)%crocus81(t)%SNOWAGE(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer(s) liquid water content (m)
            do l=1, CROCUS81_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWLIQ", &
                                         dim=l, vlevels = CROCUS81_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    CROCUS81_struc(n)%crocus81(t)%SNOWLIQ(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer(s) temperature (K)
            do l=1, CROCUS81_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWTEMP", &
                                         dim=l, vlevels = CROCUS81_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    CROCUS81_struc(n)%crocus81(t)%SNOWTEMP(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Snow layer(s) thickness (m)
            do l=1, CROCUS81_struc(n)%nsnow ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SNOWDZ", &
                                         dim=l, vlevels = CROCUS81_struc(n)%nsnow, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    CROCUS81_struc(n)%crocus81(t)%SNOWDZ(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: Soil/snow interface heat flux (W/m2)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%GRNDFLUX, &
                                     varname="GRNDFLUX", wformat=wformat)
 
            ! read: Blowing snow sublimation (kg/m2/s) NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%SNDRIFT, &
                                     varname="SNDRIFT", wformat=wformat)
 
            ! read: Richardson number (-)  NOTE: RI has not been initialized in CALL_MODEL (If not OMED initalized to undefined in the snow3L_isba.F90)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%RI_n, &
                                     varname="RI_N", wformat=wformat)
 
            ! read: Drag coefficient for momentum over snow (-)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%CDSNOW, &
                                     varname="CDSNOW", wformat=wformat)
 
            ! read: Friction velocity over snow (m/s);
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%USTARSNOW, &
                                     varname="USTARSNOW", wformat=wformat)
 
            ! read: Drag coefficient for heat over snow  (-)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%CHSNOW, &
                                     varname="CHSNOW", wformat=wformat)
 
            ! read: Snowmaking thickness (m)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%SNOWMAK_dz, &
                                     varname="SNOWMAK_DZ", wformat=wformat)
        
            ! close restart file
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in Crocus81_readrst")
#endif
            endif
            deallocate(tmptilen)
        endif    
    enddo
end subroutine Crocus81_readrst
        
