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
! !ROUTINE: AWRAL600_readrst
! \label{AWRAL600_readrst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   12/18/18: Wendy Sharples, Shugong Wang; initial implementation for LIS 7 and AWRAL600
!
! !INTERFACE:
subroutine AWRAL600_readrst()
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify                
    use AWRAL600_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for AWRAL600.  This
!  includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the AWRAL600
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    sr                         - AWRAL600 volume of water in the surface water store [mm]
!    sg                         - AWRAL600 groundwater storage in the unconfined aquifer [mm]
!    s0                         - AWRAL600 water storage in the surface soil layer for each hru [mm]
!    ss                         - AWRAL600 water content of the shallow soil store for each hru [mm]
!    sd                         - AWRAL600 water content of the deep soil store for each hru [mm]
!    mleaf                      - AWRAL600 leaf biomass [kg/m2]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart})\\
!      reads a variable from the restart file
!   \item[AWRAL600\_coldstart](\ref{AWRAL600_coldstart})\\
!      initializes the AWRAL600 state variables
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
        wformat = trim(AWRAL600_struc(n)%rformat)
        ! coldstart
        if(LIS_rc%startcode .eq. "coldstart") then  
            call AWRAL600_coldstart(LIS_rc%lsm_index)
        ! restart
        elseif(LIS_rc%startcode .eq. "restart") then
            allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            ! check the existance of restart file
            inquire(file=AWRAL600_struc(n)%rfile, exist=file_exists)
            If (.not. file_exists) then 
                write(LIS_logunit,*) "[ERR] AWRAL600 restart file ", AWRAL600_struc(n)%rfile," does not exist "
                write(LIS_logunit,*) "[ERR] Program stopping ..."
                call LIS_endrun
            endif
            write(LIS_logunit,*) "[INFO] AWRAL600 restart file used: ", AWRAL600_struc(n)%rfile
        
            ! open restart file
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn, file=AWRAL600_struc(n)%rfile, form="unformatted")
                read(ftn) nc, nr, npatch  !time, veg class, no. tiles
 
                ! check for grid space conflict
                if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
                    write(LIS_logunit,*) AWRAL600_struc(n)%rfile, "[ERR] grid space mismatch - AWRAL600 halted"
                    call LIS_endrun
                endif
            
                if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
                    write(LIS_logunit,*) "[ERR] restart tile space mismatch, halting..."
                    call LIS_endrun
                endif
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_open(path=AWRAL600_struc(n)%rfile, &
                                   mode=NF90_NOWRITE, ncid=ftn)
                call LIS_verify(status, "[ERR] Error opening file "//AWRAL600_struc(n)%rfile)
#endif            
            endif
 
            ! read: volume of water in the surface water store
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AWRAL600_struc(n)%awral600%sr, &
                                     varname="SR", wformat=wformat)
 
            ! read: groundwater storage in the unconfined aquifer
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, AWRAL600_struc(n)%awral600%sg, &
                                     varname="SG", wformat=wformat)
 
            ! read: water storage in the surface soil layer for each hru
            do l=1, AWRAL600_struc(n)%nhru
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="S0", &
                                         dim=l, vlevels = AWRAL600_struc(n)%nhru, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    AWRAL600_struc(n)%awral600(t)%s0(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: water content of the shallow soil store for each hru
            do l=1, AWRAL600_struc(n)%nhru
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SS", &
                                         dim=l, vlevels = AWRAL600_struc(n)%nhru, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    AWRAL600_struc(n)%awral600(t)%ss(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: water content of the deep soil store for each hru
            do l=1, AWRAL600_struc(n)%nhru
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SD", &
                                         dim=l, vlevels = AWRAL600_struc(n)%nhru, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    AWRAL600_struc(n)%awral600(t)%sd(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: leaf biomass
            do l=1, AWRAL600_struc(n)%nhru
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="MLEAF", &
                                         dim=l, vlevels = AWRAL600_struc(n)%nhru, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    AWRAL600_struc(n)%awral600(t)%mleaf(l) = tmptilen(t)
                enddo
            enddo
        
            ! close restart file
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in AWRAL600_readrst")
#endif            
            endif
            deallocate(tmptilen)
        endif    
    enddo
end subroutine AWRAL600_readrst
        
