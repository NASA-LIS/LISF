!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: Ac71_writerst
! \label{Ac71_writerst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   9/4/14: Shugong Wang; initial implementation for LIS 7 and Ac71
!   LB: TODO
!
! !INTERFACE:
subroutine Ac71_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use Ac71_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for Ac71.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[Ac71\_dump\_restart](\ref{Ac71_dump_restart}) \newline
!   writes the Ac71 variables into the restart file
! \end{description}
!EOP

    character(len=LIS_CONST_PATH_LEN) :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn
    integer       :: status
    
    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "Ac71 restart alarm")
    
    ! set restart file format (read from LIS configration file_
    wformat = trim(AC71_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", "AC71",&
                                             wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status,"Error in nf90_open in Ac71_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, "Error in nf90_open in Ac71_writerst")
#endif
             endif
        endif
    
        call Ac71_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in Ac71_writerst")
#endif
            endif
            write(LIS_logunit, *) "Ac71 archive restart written: ", trim(filen)
        endif
    endif
end subroutine Ac71_writerst

!BOP
!
! !ROUTINE: Ac71_dump_restart
! \label{Ac71_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  9/4/14: Shugong Wang, initial implementation for LIS 7 and Ac71
! !INTERFACE:
subroutine Ac71_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use Ac71_lsmMod

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
!  The following is the list of variables written in the Ac71
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
         "AC71", &
         dim1 = AC71_struc(n)%nsoil, &
         dimID = dimID, &
         output_format = trim(wformat))

end subroutine Ac71_dump_restart
