!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: noahmpglacier3911_writerst
! \label{noahmpglacier3911_writerst}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine noahmpglacier3911_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use noahmpglacier3911_Mod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for noahmpglacier3911.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[noahmpglacier3911\_dump\_restart](\ref{noahmpglacier3911_dump_restart}) \newline
!   writes the noahmpglacier3911 variables into the restart file
! \end{description}
!EOP

    character*100 :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn
    integer       :: status

    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "noahmpglacier3911 restart alarm")
    
    ! set restart file format (read from LIS configration file_
    wformat = trim(noahmpgl3911_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", "NOAHMPGL3911",&
                                             wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status,"Error in nf90_open in noahmpglacier3911_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, "Error in nf90_open in noahmpglacier3911_writerst")
#endif
             endif
        endif
    
        call noahmpglacier3911_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in noahmpglacier3911_writerst")
#endif
            endif
            write(LIS_logunit, *) "noahmpglacier3911 archive restart written: ", filen
        endif
    endif

end subroutine noahmpglacier3911_writerst

!BOP
!
! !ROUTINE: noahmpglacier3911_dump_restart
! \label{noahmpglacier3911_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  9/4/14: Shugong Wang, initial implementation for LIS 7 and noahmpglacier3911
! !INTERFACE:
subroutine noahmpglacier3911_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use noahmpglacier3911_Mod

    implicit none

    integer, intent(in) :: ftn
    integer, intent(in) :: n
    character(len=*), intent(in) :: wformat
!
! !DESCRIPTION:
!  This routine gathers the necessary restart variables and performs
!  the actual write statements to create the restart files.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    unit number for the restart file
!   \item[wformat]
!    restart file format (binary/netcdf)
!  \end{description}
!
!
!  The following is the list of variables written in the noahmpglacier3911
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
! The routines invoked are:
! \begin{description}
!   \item[LIS\_writeGlobalHeader\_restart](\ref{LIS_writeGlobalHeader_restart}) \newline
!      writes the global header information
!   \item[LIS\_writeHeader\_restart](\ref{LIS_writeHeader_restart}) \newline
!      writes the header information for a variable
!   \item[LIS\_closeHeader\_restart](\ref{LIS_closeHeader_restart}) \newline
!      close the header
!   \item[LIS\_writevar\_restart](\ref{LIS_writevar_restart}) \newline
!      writes a variable to the restart file
! \end{description}
! 
!EOP 
               
    integer :: l, t 
    real    :: tmptilen(LIS_rc%npatch(n, LIS_rc%glacier_index))
    integer :: dimID(10)
    integer :: albold_ID
    integer :: sneqvo_ID
    integer :: sstc_ID
    integer :: sh2o_ID
    integer :: smc_ID
    integer :: tah_ID
    integer :: eah_ID
    integer :: fwet_ID
    integer :: canliq_ID
    integer :: canice_ID
    integer :: tv_ID
    integer :: tg_ID
    integer :: qsnow_ID
    integer :: isnow_ID
    integer :: zss_ID
    integer :: snowh_ID
    integer :: sneqv_ID
    integer :: snowice_ID
    integer :: snowliq_ID
    integer :: zwt_ID
    integer :: wa_ID
    integer :: wt_ID
    integer :: wslake_ID
    integer :: lfmass_ID
    integer :: rtmass_ID
    integer :: stmass_ID
    integer :: wood_ID
    integer :: stblcp_ID
    integer :: fastcp_ID
    integer :: lai_ID
    integer :: sai_ID
    integer :: cm_ID
    integer :: ch_ID
    integer :: tauss_ID
    integer :: smcwtd_ID
    integer :: deeprech_ID
    integer :: rech_ID
    integer :: zlvl_ID
    
    ! write the header of the restart file
    call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%glacier_index, &
         "NOAHMP36", &
         dim1 = noahmpgl3911_struc(n)%nsoil+noahmpgl3911_struc(n)%nsnow,   &
         dim2 = noahmpgl3911_struc(n)%nsoil, &
         dim3 = noahmpgl3911_struc(n)%nsnow, & 
         dimID = dimID, &
         output_format = trim(wformat))
    
    ! write the header for state variable albold
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, albold_ID, "ALBOLD", &
         "snow albedo at last time step", &
         "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable sneqvo
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sneqvo_ID, "SNEQVO", &
         "snow mass at the last time step", &
         "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable sstc
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sstc_ID, "SSTC", &
         "snow/soil temperature", &
         "K", vlevels=noahmpgl3911_struc(n)%nsoil+noahmpgl3911_struc(n)%nsnow, &
         valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable sh2o
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sh2o_ID, "SH2O", &
         "volumetric liquid soil moisture", &
         "m^3 m-3", vlevels=noahmpgl3911_struc(n)%nsoil, valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim2") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable smc
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, smc_ID, "SMC", &
         "volumetric soil moisture, ice + liquid", &
         "m^3 m-3", vlevels=noahmpgl3911_struc(n)%nsoil , valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim2") ! TODO: replace "xxxx" with correct "dimx"
 
    call LIS_writeHeader_restart(ftn, n, dimID, tg_ID, "TG", &
         "ground temperature (skin temperature)", &
         "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable qsnow
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qsnow_ID, "QSNOW", &
         "snowfall on the ground", &
         "mm s-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable isnow
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, isnow_ID, "ISNOW", &
         "actual number of snow layers", &
         "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable zss
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, zss_ID, "ZSS", &
         "snow/soil layer-bottom depth from snow surface", &
         "m", vlevels=noahmpgl3911_struc(n)%nsoil+noahmpgl3911_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim1") ! TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable snowh
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, snowh_ID, "SNOWH", &
         "snow height", &
         "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable sneqv
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sneqv_ID, "SNEQV", &
         "snow water equivalent", &
         "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable snowice
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, snowice_ID, "SNOWICE", &
         "snow-layer ice", &
         "mm", vlevels=noahmpgl3911_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim3") ! TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable snowliq
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, snowliq_ID, "SNOWLIQ", &
         "snow-layer liquid water", &
         "mm", vlevels=noahmpgl3911_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim3") ! TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable cm
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, cm_ID, "CM", &
         "momentum drag coefficient", &
         "s m-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable ch
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, ch_ID, "CH", &
         "sensible heat exchange coefficient", &
         "s m-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable tauss
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, tauss_ID, "TAUSS", &
         "snow aging term", &
         "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable zlvl
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, zlvl_ID, "ZLVL", &
         "reference height for air temperature and humidity", &
         "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! close header of restart file
    call LIS_closeHeader_restart(ftn, n, LIS_rc%glacier_index, dimID, noahmpgl3911_struc(n)%rstInterval)
    
    ! write state variables into restart file
    ! snow albedo at last time step
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%albold, &
         varid=albold_ID, dim=1, wformat=wformat)
    
    ! snow mass at the last time step
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%sneqvo, &
         varid=sneqvo_ID, dim=1, wformat=wformat)

    ! snow/soil temperature
    do l=1, noahmpgl3911_struc(n)%nsoil+noahmpgl3911_struc(n)%nsnow   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
          tmptilen(t) = noahmpgl3911_struc(n)%noahmpgl(t)%sstc(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, &
            varid=sstc_ID, dim=l, wformat=wformat)
    enddo
    ! volumetric liquid soil moisture
    do l=1, noahmpgl3911_struc(n)%nsoil   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
          tmptilen(t) = noahmpgl3911_struc(n)%noahmpgl(t)%sh2o(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, &
            varid=sh2o_ID, dim=l, wformat=wformat)
    enddo
    ! volumetric soil moisture, ice + liquid
    do l=1, noahmpgl3911_struc(n)%nsoil   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
          tmptilen(t) = noahmpgl3911_struc(n)%noahmpgl(t)%smc(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, &
            varid=smc_ID, dim=l, wformat=wformat)
    enddo

    ! ground temperature (skin temperature)
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%tg, &
         varid=tg_ID, dim=1, wformat=wformat)

    ! snowfall on the ground
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%qsnow, &
         varid=qsnow_ID, dim=1, wformat=wformat)

    ! actual number of snow layers
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%isnow, &
         varid=isnow_ID, dim=1, wformat=wformat)

    ! snow/soil layer-bottom depth from snow surface
    do l=1, noahmpgl3911_struc(n)%nsoil+noahmpgl3911_struc(n)%nsnow   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
          tmptilen(t) = noahmpgl3911_struc(n)%noahmpgl(t)%zss(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, &
            varid=zss_ID, dim=l, wformat=wformat)
    enddo

    ! snow height
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%snowh, &
         varid=snowh_ID, dim=1, wformat=wformat)

    ! snow water equivalent
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%sneqv, &
         varid=sneqv_ID, dim=1, wformat=wformat)

    ! snow-layer ice
    do l=1, noahmpgl3911_struc(n)%nsnow   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
          tmptilen(t) = noahmpgl3911_struc(n)%noahmpgl(t)%snowice(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, &
            varid=snowice_ID, dim=l, wformat=wformat)
    enddo
    ! snow-layer liquid water
    do l=1, noahmpgl3911_struc(n)%nsnow   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%glacier_index)
          tmptilen(t) = noahmpgl3911_struc(n)%noahmpgl(t)%snowliq(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, tmptilen, &
            varid=snowliq_ID, dim=l, wformat=wformat)
    enddo

    ! momentum drag coefficient
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%cm, &
         varid=cm_ID, dim=1, wformat=wformat)

    ! sensible heat exchange coefficient
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%ch, &
         varid=ch_ID, dim=1, wformat=wformat)

    ! snow aging term
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%tauss, &
         varid=tauss_ID, dim=1, wformat=wformat)

    ! reference height for air temperature and humidity 
    call LIS_writevar_restart(ftn, n, LIS_rc%glacier_index, noahmpgl3911_struc(n)%noahmpgl%zlvl, &
         varid=zlvl_ID, dim=1, wformat=wformat)

end subroutine noahmpglacier3911_dump_restart
