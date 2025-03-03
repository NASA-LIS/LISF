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
! !ROUTINE: noahmpglacier3911_readrst
! \label{noahmpglacier3911_readrst}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine noahmpglacier3911_readrst()
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_historyMod, only : LIS_readvar_restart
    use LIS_logMod, only     : LIS_logunit, LIS_endrun, &
                               LIS_getNextUnitNumber,   &
                               LIS_releaseUnitNumber,   &
                               LIS_verify                
    use noahmpglacier3911_Mod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

!
! !DESCRIPTION:
!  This program reads restart files for noahmpglacier3911.  This
!  includes all relevant water/energy storages and tile information.
!  The following is the list of variables specified in the noahmpglacier3911
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    albold                     - noahmpglacier3911 snow albedo at last time step [-]
!    sneqvo                     - noahmpglacier3911 snow mass at the last time step [mm]
!    sstc                       - noahmpglacier3911 snow/soil temperature [K]
!    sh2o                       - noahmpglacier3911 volumetric liquid soil moisture [m^3 m-3]
!    smc                        - noahmpglacier3911 volumetric soil moisture, ice + liquid [m^3 m-3]
!    tah                        - noahmpglacier3911 canopy air temperature [K]
!    eah                        - noahmpglacier3911 canopy air vapor pressure [Pa]
!    fwet                       - noahmpglacier3911 wetted or snowed fraction of canopy [-]
!    canliq                     - noahmpglacier3911 intercepted liquid water [mm]
!    canice                     - noahmpglacier3911 intercepted ice mass [mm]
!    tv                         - noahmpglacier3911 vegetation temperature [K]
!    tg                         - noahmpglacier3911 ground temperature (skin temperature) [K]
!    qsnow                      - noahmpglacier3911 snowfall on the ground [mm s-1]
!    isnow                      - noahmpglacier3911 actual number of snow layers [-]
!    zss                        - noahmpglacier3911 snow/soil layer-bottom depth from snow surface [m]
!    snowh                      - noahmpglacier3911 snow height [m]
!    sneqv                      - noahmpglacier3911 snow water equivalent [mm]
!    snowice                    - noahmpglacier3911 snow-layer ice [mm]
!    snowliq                    - noahmpglacier3911 snow-layer liquid water [mm]
!    zwt                        - noahmpglacier3911 depth to water table [m]
!    wa                         - noahmpglacier3911 water storage in aquifer [mm]
!    wt                         - noahmpglacier3911 water in aquifer and saturated soil [mm]
!    wslake                     - noahmpglacier3911 lake water storage [mm]
!    lfmass                     - noahmpglacier3911 leaf mass [g/m2]
!    rtmass                     - noahmpglacier3911 mass of fine roots [g/m2]
!    stmass                     - noahmpglacier3911 stem mass [g/m2]
!    wood                       - noahmpglacier3911 mass of wood including woody roots [g/m2]
!    stblcp                     - noahmpglacier3911 stable carbon in deep soil [g/m2]
!    fastcp                     - noahmpglacier3911 short-lived carbon in shallow soil [g/m2]
!    lai                        - noahmpglacier3911 leaf area index [-]
!    sai                        - noahmpglacier3911 stem area index [-]
!    cm                         - noahmpglacier3911 momentum drag coefficient [s/m]
!    ch                         - noahmpglacier3911 sensible heat exchange coefficient [s/m]
!    tauss                      - noahmpglacier3911 snow aging term [-]
!    smcwtd                     - noahmpglacier3911 soil water content between bottom of the soil and water table [m^3 m-3]
!    deeprech                   - noahmpglacier3911 recharge to or from the water table when deep [m]
!    rech                       - noahmpglacier3911 recharge to or from the water table when shallow [m]
!  \end{verbatim}
!
!  The routines invoked are:
! \begin{description}
!   \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!      reads a variable from the restart file
!   \item[noahmpglacier3911\_coldstart](\ref{noahmpglacier3911_coldstart}) \newline
!      initializes the noahmpglacier3911 state variables
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
       wformat = trim(noahmpgl3911_struc(n)%rformat)
       ! coldstart
       if(LIS_rc%startcode .eq. "coldstart") then  
          call noahmpglacier3911_coldstart(LIS_rc%glacier_index)
          ! restart
       elseif(LIS_rc%startcode .eq. "restart") then
          allocate(tmptilen(LIS_rc%npatch(n, LIS_rc%glacier_index)))
          ! check the existance of restart file
          inquire(file=noahmpgl3911_struc(n)%rfile, exist=file_exists)
          If (.not. file_exists) then 
             write(LIS_logunit,*) "noahmpglacier3911 restart file ", noahmpgl3911_struc(n)%rfile," does not exist "
             write(LIS_logunit,*) "Program stopping ..."
             call LIS_endrun
          endif
          write(LIS_logunit,*) "noahmpglacier3911 restart file used: ", noahmpgl3911_struc(n)%rfile

          ! open restart file
          if(wformat .eq. "binary") then
             ftn = LIS_getNextUnitNumber()
             open(ftn, file=noahmpgl3911_struc(n)%rfile, form="unformatted")
             read(ftn) nc, nr, npatch  !time, veg class, no. tiles

             ! check for grid space conflict
             if((nc .ne. LIS_rc%gnc(n)) .or. (nr .ne. LIS_rc%gnr(n))) then
                write(LIS_logunit,*) noahmpgl3911_struc(n)%rfile, "grid space mismatch - noahmpglacier3911 halted"
                call LIS_endrun
             endif

             if(npatch .ne. LIS_rc%glbnpatch_red(n, LIS_rc%glacier_index)) then
                write(LIS_logunit,*) "restart tile space mismatch, halting..."
                call LIS_endrun
             endif
          elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
             status = nf90_open(path=noahmpgl3911_struc(n)%rfile, &
                  mode=NF90_NOWRITE, ncid=ftn)
             call LIS_verify(status, "Error opening file "//noahmpgl3911_struc(n)%rfile)
#endif
          endif

          ! read: snow albedo at last time step
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%albold, &
               varname="ALBOLD", wformat=wformat)

          ! read: snow mass at the last time step
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%sneqvo, &
               varname="SNEQVO", wformat=wformat)

          ! read: snow/soil temperature
          do l=1, noahmpgl3911_struc(n)%nsoil + noahmpgl3911_struc(n)%nsnow
             call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, varname="SSTC", &
                  dim=l, vlevels = noahmpgl3911_struc(n)%nsoil + noahmpgl3911_struc(n)%nsnow,&
                  wformat=wformat)
             do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
                noahmpgl3911_struc(n)%noahmpgl(t)%sstc(l) = tmptilen(t)
             enddo
          enddo

          ! read: volumetric liquid soil moisture
          do l=1, noahmpgl3911_struc(n)%nsoil ! TODO: check loop
             call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, varname="SH2O", &
                  dim=l, vlevels = noahmpgl3911_struc(n)%nsoil, wformat=wformat)
             do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
                noahmpgl3911_struc(n)%noahmpgl(t)%sh2o(l) = tmptilen(t)
             enddo
          enddo

          ! read: volumetric soil moisture, ice + liquid
          do l=1, noahmpgl3911_struc(n)%nsoil ! TODO: check loop
             call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, varname="SMC", &
                  dim=l, vlevels = noahmpgl3911_struc(n)%nsoil, wformat=wformat)
             do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
                noahmpgl3911_struc(n)%noahmpgl(t)%smc(l) = tmptilen(t)
             enddo
          enddo

          ! read: ground temperature (skin temperature)
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%tg, &
               varname="TG", wformat=wformat)

          ! read: snowfall on the ground
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%qsnow, &
               varname="QSNOW", wformat=wformat)

          ! read: actual number of snow layers
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%isnow, &
               varname="ISNOW", wformat=wformat)

          ! read: snow/soil layer-bottom depth from snow surface
          do l=1, noahmpgl3911_struc(n)%nsoil + noahmpgl3911_struc(n)%nsnow
             call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, varname="ZSS", &
                  dim=l, vlevels = noahmpgl3911_struc(n)%nsoil + noahmpgl3911_struc(n)%nsnow, &
                  wformat=wformat)
             do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
                noahmpgl3911_struc(n)%noahmpgl(t)%zss(l) = tmptilen(t)
             enddo
          enddo

          ! read: snow height
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%snowh, &
               varname="SNOWH", wformat=wformat)

          ! read: snow water equivalent
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%sneqv, &
               varname="SNEQV", wformat=wformat)

          ! read: snow-layer ice
          do l=1, noahmpgl3911_struc(n)%nsnow ! TODO: check loop
             call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, varname="SNOWICE", &
                  dim=l, vlevels = noahmpgl3911_struc(n)%nsnow, wformat=wformat)
             do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
                noahmpgl3911_struc(n)%noahmpgl(t)%snowice(l) = tmptilen(t)
             enddo
          enddo

          ! read: snow-layer liquid water
          do l=1, noahmpgl3911_struc(n)%nsnow ! TODO: check loop
             call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, varname="SNOWLIQ", &
                  dim=l, vlevels = noahmpgl3911_struc(n)%nsnow, wformat=wformat)
             do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
                noahmpgl3911_struc(n)%noahmpgl(t)%snowliq(l) = tmptilen(t)
             enddo
          enddo

          ! read: momentum drag coefficient
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%cm, &
               varname="CM", wformat=wformat)

          ! read: sensible heat exchange coefficient
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%ch, &
               varname="CH", wformat=wformat)

          ! read: snow aging term
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%tauss, &
               varname="TAUSS", wformat=wformat)

          ! read: reference height for air temperature and humidity 
          call LIS_readvar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%zlvl, &
               varname="ZLVL", wformat=wformat)
          ! close restart file
          if(wformat .eq. "binary") then
             call LIS_releaseUnitNumber(ftn)
          elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
             status = nf90_close(ftn)
             call LIS_verify(status, "Error in nf90_close in noahmpglacier3911_readrst")
#endif
          endif
          deallocate(tmptilen)

       endif

    enddo
    
  end subroutine noahmpglacier3911_readrst
        
