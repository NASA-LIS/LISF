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
! !ROUTINE: RUC37_readrst
! \label{RUC37_readrst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   1/15/15: Shugong Wang; initial implementation for LIS 7 and RUC37
!
! !INTERFACE:
subroutine RUC37_readrst()
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify                
    use RUC37_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for RUC37.  This
!  includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the RUC37
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    emiss                      - RUC37 surface emissivity (0.0 - 1.0). [-]
!    ch                         - RUC37 exchange coefficient for head and moisture (m s-1). [s/m]
!    cm                         - RUC37 exchange coefficient for momentum (m s-1). [s/m]
!    sneqv                      - RUC37 water equivalent of accumulated snow depth (m). [m]
!    snowh                      - RUC37 physical snow depth (m). [m]
!    snowc                      - RUC37 fractional snow cover ( fraction [0.0-1.0] ) [-]
!    canwat                     - RUC37 canopy moisture content (kg m-2) [kg m-2]
!    alb                        - RUC37 surface albedo including possible snow-cover effect.  this is set in lsmruc, [-]
!    smc                        - RUC37 total soil moisture content (m3 m-3) [m^3 m-3]
!    sho                        - RUC37 liquid soil moisture content (m3 m-3) [m^3 m-3]
!    stc                        - RUC37 soil temperature (k) [K]
!    smfr                       - RUC37 soil ice (m3 m-3) [m^3 m-3]
!    keepfr                     - RUC37 frozen soil flag
!    tskin                      - RUC37 skin temperature (k) [K]
!    qvg                        - RUC37 mixing ratio at the surface ( kg kg{-1} ) [kg kg-1]
!    qsfc                       - RUC37 specific humidity at the surface ( kg kg{-1} ) [kg kg-1]
!    qcg                        - RUC37 cloud water mixing ratio at the surface ( kg kg{-1} ) [kg kg-1]
!    qsg                        - RUC37 surface water vapor mixing ratio at satration (kg kg-1) [kg/kg]
!    snt75cm                    - RUC37 snow temperature at 7.5 cm depth (k) [K]
!    tsnav                      - RUC37 average snow temperature in k [K]
!    soilm                      - RUC37 total soil column moisture content, frozen and unfrozen ( m ) [m]
!    smroot                     - RUC37 available soil moisture in the root zone ( fraction [smcwlt-smcmax] [m^3 m-3]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!      reads a variable from the restart file
!   \item[RUC37\_coldstart](\ref{RUC37_coldstart}) \newline
!      initializes the RUC37 state variables
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
        wformat = trim(RUC37_struc(n)%rformat)
        ! coldstart
        if(LIS_rc%startcode .eq. "coldstart") then  
            call RUC37_coldstart(LIS_rc%lsm_index)
        ! restart
        elseif(LIS_rc%startcode .eq. "restart") then
            allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            ! check the existance of restart file
            inquire(file=RUC37_struc(n)%rfile, exist=file_exists)
            If (.not. file_exists) then 
                write(LIS_logunit,*) "RUC37 restart file ", RUC37_struc(n)%rfile," does not exist "
                write(LIS_logunit,*) "Program stopping ..."
                call LIS_endrun
            endif
            write(LIS_logunit,*) "RUC37 restart file used: ", RUC37_struc(n)%rfile
        
            ! open restart file
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn, file=RUC37_struc(n)%rfile, form="unformatted")
                read(ftn) nc, nr, npatch  !time, veg class, no. tiles
 
                ! check for grid space conflict
                if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
                    write(LIS_logunit,*) RUC37_struc(n)%rfile, "grid space mismatch - RUC37 halted"
                    call LIS_endrun
                endif
            
                if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%lsm_index)) then
                    write(LIS_logunit,*) "restart tile space mismatch, halting..."
                    call LIS_endrun
                endif
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 ||  defined USE_NETCDF4)
                status = nf90_open(path=RUC37_struc(n)%rfile, &
                                   mode=NF90_NOWRITE, ncid=ftn)
                call LIS_verify(status, "Error opening file "//RUC37_struc(n)%rfile)
#endif
            endif
 
            ! read: surface emissivity (0.0 - 1.0).
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%emiss, &
                                     varname="EMISS", wformat=wformat)
 
            ! read: exchange coefficient for head and moisture (m s-1).
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%ch, &
                                     varname="CH", wformat=wformat)
 
            ! read: exchange coefficient for momentum (m s-1).
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%cm, &
                                     varname="CM", wformat=wformat)
 
            ! read: water equivalent of accumulated snow depth (m).
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%sneqv, &
                                     varname="SNEQV", wformat=wformat)
 
            ! read: physical snow depth (m).
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%snowh, &
                                     varname="SNOWH", wformat=wformat)
 
            ! read: fractional snow cover ( fraction [0.0-1.0] )
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%snowc, &
                                     varname="SNOWC", wformat=wformat)
 
            ! read: canopy moisture content (kg m-2)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%canwat, &
                                     varname="CANWAT", wformat=wformat)
 
            ! read: surface albedo including possible snow-cover effect.  this is set in lsmruc,
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%alb, &
                                     varname="ALB", wformat=wformat)
 
            ! read: total soil moisture content (m3 m-3)
            do l=1, RUC37_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SMC", &
                                         dim=l, vlevels = RUC37_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    RUC37_struc(n)%ruc37(t)%smc(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: liquid soil moisture content (m3 m-3)
            do l=1, RUC37_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SHO", &
                                         dim=l, vlevels = RUC37_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    RUC37_struc(n)%ruc37(t)%sho(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: soil temperature (k)
            do l=1, RUC37_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="STC", &
                                         dim=l, vlevels = RUC37_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    RUC37_struc(n)%ruc37(t)%stc(l) = tmptilen(t)
                enddo
            enddo
 
            ! read: soil ice content (m3 m-3)
            do l=1, RUC37_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="SMFR", &
                                         dim=l, vlevels = RUC37_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    RUC37_struc(n)%ruc37(t)%smfr(l) = tmptilen(t)
                enddo
            enddo

            ! read: frozen soil flag
            do l=1, RUC37_struc(n)%nsoil ! TODO: check loop
                call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="KEEPFR", &
                                         dim=l, vlevels = RUC37_struc(n)%nsoil, wformat=wformat)
                do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                    RUC37_struc(n)%ruc37(t)%keepfr(l) = tmptilen(t)
                enddo
            enddo

            ! read: skin temperature (k)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%tskin, &
                                     varname="TSKIN", wformat=wformat)
 
            ! read: mixing ratio at the surface ( kg kg{-1} )
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%qvg, &
                                     varname="QVG", wformat=wformat)

            ! read: specific humidity at the surface ( kg kg{-1} )
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%qsfc, &
                                     varname="QSFC", wformat=wformat)

            ! read: cloud water mixing ratio at the surface ( kg kg{-1} )
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%qcg, &
                                     varname="QCG", wformat=wformat)
 
            ! read: surface water vapor mixing ratio at satration (kg kg-1)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%qsg, &
                                     varname="QSG", wformat=wformat)
 
            ! read: snow temperature at 7.5 cm depth (k)
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%snt75cm, &
                                     varname="SNT75CM", wformat=wformat)
 
            ! read: average snow temperature in k
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%tsnav, &
                                     varname="TSNAV", wformat=wformat)
 
            ! read: total soil column moisture content, frozen and unfrozen ( m )
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%soilm, &
                                     varname="SOILM", wformat=wformat)
 
            ! read: available soil moisture in the root zone ( fraction [smcwlt-smcmax]
            call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%smroot, &
                                     varname="SMROOT", wformat=wformat)
        
            ! close restart file
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in RUC37_readrst")
#endif
            endif
            deallocate(tmptilen)
        endif    
    enddo
end subroutine RUC37_readrst
        
