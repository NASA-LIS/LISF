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
! !ROUTINE: NoahMP36_writerst
! \label{NoahMP36_writerst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   9/4/14: Shugong Wang; initial implementation for LIS 7 and NoahMP36
!
! !INTERFACE:
subroutine NoahMP36_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use NoahMP36_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for NoahMP36.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[NoahMP36\_dump\_restart](\ref{NoahMP36_dump_restart}) \newline
!   writes the NoahMP36 variables into the restart file
! \end{description}
!EOP

    character(len=LIS_CONST_PATH_LEN) :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn
    integer       :: status
    
    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "NoahMP36 restart alarm")
    
    ! set restart file format (read from LIS configration file_
    wformat = trim(NOAHMP36_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", "NOAHMP36",&
                                             wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status,"Error in nf90_open in NoahMP36_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, "Error in nf90_open in NoahMP36_writerst")
#endif
             endif
        endif
    
        call NoahMP36_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in NoahMP36_writerst")
#endif
            endif
            write(LIS_logunit, *) "NoahMP36 archive restart written: ", trim(filen)
        endif
    endif
end subroutine NoahMP36_writerst

!BOP
!
! !ROUTINE: NoahMP36_dump_restart
! \label{NoahMP36_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  9/4/14: Shugong Wang, initial implementation for LIS 7 and NoahMP36
! !INTERFACE:
subroutine NoahMP36_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use NoahMP36_lsmMod

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
!  The following is the list of variables written in the NoahMP36
!  restart file:
!  \begin{verbatim}
!    nc, nr, ntiles             - grid and tile space dimensions
!    albold                     - NoahMP36 snow albedo at last time step [-]
!    sneqvo                     - NoahMP36 snow mass at the last time step [mm]
!    sstc                       - NoahMP36 snow/soil temperature [K]
!    sh2o                       - NoahMP36 volumetric liquid soil moisture [m^3 m-3]
!    smc                        - NoahMP36 volumetric soil moisture, ice + liquid [m^3 m-3]
!    tah                        - NoahMP36 canopy air temperature [K]
!    eah                        - NoahMP36 canopy air vapor pressure [Pa]
!    fwet                       - NoahMP36 wetted or snowed fraction of canopy [-]
!    canliq                     - NoahMP36 intercepted liquid water [mm]
!    canice                     - NoahMP36 intercepted ice mass [mm]
!    tv                         - NoahMP36 vegetation temperature [K]
!    tg                         - NoahMP36 ground temperature (skin temperature) [K]
!    qsnow                      - NoahMP36 snowfall on the ground [mm s-1]
!    isnow                      - NoahMP36 actual number of snow layers [-]
!    zss                        - NoahMP36 snow/soil layer-bottom depth from snow surface [m]
!    snowh                      - NoahMP36 snow height [m]
!    sneqv                      - NoahMP36 snow water equivalent [mm]
!    snowice                    - NoahMP36 snow-layer ice [mm]
!    snowliq                    - NoahMP36 snow-layer liquid water [mm]
!    zwt                        - NoahMP36 depth to water table [m]
!    wa                         - NoahMP36 water storage in aquifer [mm]
!    wt                         - NoahMP36 water in aquifer and saturated soil [mm]
!    wslake                     - NoahMP36 lake water storage [mm]
!    lfmass                     - NoahMP36 leaf mass [g/m2]
!    rtmass                     - NoahMP36 mass of fine roots [g/m2]
!    stmass                     - NoahMP36 stem mass [g/m2]
!    wood                       - NoahMP36 mass of wood including woody roots [g/m2]
!    stblcp                     - NoahMP36 stable carbon in deep soil [g/m2]
!    fastcp                     - NoahMP36 short-lived carbon in shallow soil [g/m2]
!    lai                        - NoahMP36 leaf area index [-]
!    sai                        - NoahMP36 stem area index [-]
!    cm                         - NoahMP36 momentum drag coefficient [s/m]
!    ch                         - NoahMP36 sensible heat exchange coefficient [s/m]
!    tauss                      - NoahMP36 snow aging term [-]
!    smcwtd                     - NoahMP36 soil water content between bottom of the soil and water table [m^3 m-3]
!    deeprech                   - NoahMP36 recharge to or from the water table when deep [m]
!    rech                       - NoahMP36 recharge to or from the water table when shallow [m]
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
    real    :: tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index))
    integer :: dimID(11)
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
    call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, &
         "NOAHMP36", &
         dim1 = NOAHMP36_struc(n)%nsoil+NOAHMP36_struc(n)%nsnow,   &
         dim2 = NOAHMP36_struc(n)%nsoil, &
         dim3 = NOAHMP36_struc(n)%nsnow, & 
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
         "K", vlevels=NOAHMP36_struc(n)%nsoil+NOAHMP36_struc(n)%nsnow, &
         valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable sh2o
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sh2o_ID, "SH2O", &
         "volumetric liquid soil moisture", &
         "m^3 m-3", vlevels=NOAHMP36_struc(n)%nsoil, valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim2") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable smc
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, smc_ID, "SMC", &
         "volumetric soil moisture, ice + liquid", &
         "m^3 m-3", vlevels=NOAHMP36_struc(n)%nsoil , valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim2") ! TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable tah
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, tah_ID, "TAH", &
         "canopy air temperature", &
         "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable eah
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, eah_ID, "EAH", &
         "canopy air vapor pressure", &
         "Pa", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable fwet
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, fwet_ID, "FWET", &
         "wetted or snowed fraction of canopy", &
         "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable canliq
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, canliq_ID, "CANLIQ", &
         "intercepted liquid water", &
         "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable canice
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, canice_ID, "CANICE", &
         "intercepted ice mass", &
         "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable tv
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, tv_ID, "TV", &
         "vegetation temperature", &
         "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable tg
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
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
         "m", vlevels=NOAHMP36_struc(n)%nsoil+NOAHMP36_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
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
         "mm", vlevels=NOAHMP36_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim3") ! TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable snowliq
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, snowliq_ID, "SNOWLIQ", &
         "snow-layer liquid water", &
         "mm", vlevels=NOAHMP36_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
         var_flag = "dim3") ! TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable zwt
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, zwt_ID, "ZWT", &
         "depth to water table", &
         "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable wa
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, wa_ID, "WA", &
         "water storage in aquifer", &
         "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable wt
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, wt_ID, "WT", &
         "water in aquifer and saturated soil", &
         "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable wslake
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, wslake_ID, "WSLAKE", &
         "lake water storage", &
         "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable lfmass
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, lfmass_ID, "LFMASS", &
         "leaf mass", &
         "g m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable rtmass
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, rtmass_ID, "RTMASS", &
         "mass of fine roots", &
         "g m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable stmass
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, stmass_ID, "STMASS", &
         "stem mass", &
         "g m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable wood
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, wood_ID, "WOOD", &
         "mass of wood including woody roots", &
         "g m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable stblcp
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, stblcp_ID, "STBLCP", &
         "stable carbon in deep soil", &
         "g m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable fastcp
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, fastcp_ID, "FASTCP", &
         "short-lived carbon in shallow soil", &
         "g m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable lai
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, lai_ID, "LAI", &
         "leaf area index", &
         "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable sai
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sai_ID, "SAI", &
         "stem area index", &
         "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
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
    ! write the header for state variable smcwtd
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, smcwtd_ID, "SMCWTD", &
         "soil water content between bottom of the soil and water table", &
         "m^3 m-3", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable deeprech
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, deeprech_ID, "DEEPRECH", &
         "recharge to or from the water table when deep", &
         "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable rech
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, rech_ID, "RECH", &
         "recharge to or from the water table when shallow", &
         "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable zlvl
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, zlvl_ID, "ZLVL", &
         "reference height for air temperature and humidity", &
         "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! close header of restart file
    call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, NOAHMP36_struc(n)%rstInterval)
    
    ! write state variables into restart file
    ! snow albedo at last time step
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%albold, &
         varid=albold_ID, dim=1, wformat=wformat)
    
    ! snow mass at the last time step
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%sneqvo, &
         varid=sneqvo_ID, dim=1, wformat=wformat)

    ! snow/soil temperature
    do l=1, NOAHMP36_struc(n)%nsoil+NOAHMP36_struc(n)%nsnow   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          tmptilen(t) = NOAHMP36_struc(n)%noahmp36(t)%sstc(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
            varid=sstc_ID, dim=l, wformat=wformat)
    enddo
    ! volumetric liquid soil moisture
    do l=1, NOAHMP36_struc(n)%nsoil   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          tmptilen(t) = NOAHMP36_struc(n)%noahmp36(t)%sh2o(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
            varid=sh2o_ID, dim=l, wformat=wformat)
    enddo
    ! volumetric soil moisture, ice + liquid
    do l=1, NOAHMP36_struc(n)%nsoil   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          tmptilen(t) = NOAHMP36_struc(n)%noahmp36(t)%smc(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
            varid=smc_ID, dim=l, wformat=wformat)
    enddo
    ! canopy air temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%tah, &
         varid=tah_ID, dim=1, wformat=wformat)

    ! canopy air vapor pressure
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%eah, &
         varid=eah_ID, dim=1, wformat=wformat)

    ! wetted or snowed fraction of canopy
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%fwet, &
         varid=fwet_ID, dim=1, wformat=wformat)

    ! intercepted liquid water
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%canliq, &
         varid=canliq_ID, dim=1, wformat=wformat)

    ! intercepted ice mass
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%canice, &
         varid=canice_ID, dim=1, wformat=wformat)

    ! vegetation temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%tv, &
         varid=tv_ID, dim=1, wformat=wformat)

    ! ground temperature (skin temperature)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%tg, &
         varid=tg_ID, dim=1, wformat=wformat)

    ! snowfall on the ground
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%qsnow, &
         varid=qsnow_ID, dim=1, wformat=wformat)

    ! actual number of snow layers
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%isnow, &
         varid=isnow_ID, dim=1, wformat=wformat)

    ! snow/soil layer-bottom depth from snow surface
    do l=1, NOAHMP36_struc(n)%nsoil+NOAHMP36_struc(n)%nsnow   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          tmptilen(t) = NOAHMP36_struc(n)%noahmp36(t)%zss(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
            varid=zss_ID, dim=l, wformat=wformat)
    enddo

    ! snow height
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%snowh, &
         varid=snowh_ID, dim=1, wformat=wformat)

    ! snow water equivalent
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%sneqv, &
         varid=sneqv_ID, dim=1, wformat=wformat)

    ! snow-layer ice
    do l=1, NOAHMP36_struc(n)%nsnow   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          tmptilen(t) = NOAHMP36_struc(n)%noahmp36(t)%snowice(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
            varid=snowice_ID, dim=l, wformat=wformat)
    enddo
    ! snow-layer liquid water
    do l=1, NOAHMP36_struc(n)%nsnow   ! TODO: check loop
       tmptilen = 0
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          tmptilen(t) = NOAHMP36_struc(n)%noahmp36(t)%snowliq(l)
       enddo
       call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
            varid=snowliq_ID, dim=l, wformat=wformat)
    enddo
    ! depth to water table
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%zwt, &
         varid=zwt_ID, dim=1, wformat=wformat)

    ! water storage in aquifer
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%wa, &
         varid=wa_ID, dim=1, wformat=wformat)

    ! water in aquifer and saturated soil
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%wt, &
         varid=wt_ID, dim=1, wformat=wformat)

    ! lake water storage
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%wslake, &
         varid=wslake_ID, dim=1, wformat=wformat)

    ! leaf mass
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%lfmass, &
         varid=lfmass_ID, dim=1, wformat=wformat)

    ! mass of fine roots
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%rtmass, &
         varid=rtmass_ID, dim=1, wformat=wformat)

    ! stem mass
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%stmass, &
         varid=stmass_ID, dim=1, wformat=wformat)

    ! mass of wood including woody roots
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%wood, &
         varid=wood_ID, dim=1, wformat=wformat)

    ! stable carbon in deep soil
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%stblcp, &
         varid=stblcp_ID, dim=1, wformat=wformat)

    ! short-lived carbon in shallow soil
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%fastcp, &
         varid=fastcp_ID, dim=1, wformat=wformat)

    ! leaf area index
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%lai, &
         varid=lai_ID, dim=1, wformat=wformat)

    ! stem area index
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%sai, &
         varid=sai_ID, dim=1, wformat=wformat)

    ! momentum drag coefficient
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%cm, &
         varid=cm_ID, dim=1, wformat=wformat)

    ! sensible heat exchange coefficient
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%ch, &
         varid=ch_ID, dim=1, wformat=wformat)

    ! snow aging term
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%tauss, &
         varid=tauss_ID, dim=1, wformat=wformat)

    ! soil water content between bottom of the soil and water table
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%smcwtd, &
         varid=smcwtd_ID, dim=1, wformat=wformat)

    ! recharge to or from the water table when deep
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%deeprech, &
         varid=deeprech_ID, dim=1, wformat=wformat)

    ! recharge to or from the water table when shallow
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%rech, &
         varid=rech_ID, dim=1, wformat=wformat)

    ! reference height for air temperature and humidity 
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NOAHMP36_struc(n)%noahmp36%zlvl, &
         varid=zlvl_ID, dim=1, wformat=wformat)

end subroutine NoahMP36_dump_restart
