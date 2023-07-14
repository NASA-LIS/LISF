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
! !ROUTINE: NoahMPnew_writerst
! \label{NoahMPnew_writerst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/25/18: Shugong Wang, Zhuo Wang; initial implementation for LIS 7 and NoahMP401
!  May 2023: Cenlin He, modified for refactored NoahMP v5 and later

! !INTERFACE:
subroutine NoahMPnew_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use NoahMPnew_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for Noah-MP-4.0.1 LSM.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory})\\
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename})\\
!  generates a timestamped restart filename
! \item[NoahMPnew\_dump\_restart](\ref{NoahMPnew_dump_restart})\\
!   writes the NoahMPnew variables into the restart file
! \end{description}
!EOP

    character(len=LIS_CONST_PATH_LEN) :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn
    integer       :: status
    
    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "NoahMPnew restart alarm")
    
    ! set restart file format (read from LIS configration file_
    wformat = trim(NoahMPnew_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", &
                                            "NOAHMPnew",wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status, &
                     "Error in nf90_open in NoahMPnew_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, &
                     "Error in nf90_open in NoahMPnew_writerst")
#endif
             endif
        endif
    
        call NoahMPnew_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, &
                     "Error in nf90_close in NoahMPnew_writerst")
#endif
            endif
            write(LIS_logunit, *)&
                 "[INFO] Noah-MP.New archive restart written: ",trim(filen)
        endif
    endif
end subroutine NoahMPnew_writerst

!BOP
!
! !ROUTINE: NoahMPnew_dump_restart
! \label{NoahMPnew_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  10/25/18: Shugong Wang, Zhuo Wang, initial implementation for LIS 7 and NoahMP401
!  May 2023: Cenlin He, modified for refactored NoahMP v5 and later
! !INTERFACE:
subroutine NoahMPnew_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use NoahMPnew_lsmMod

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
!  The following is the list of variables written in the NoahMP401
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
! The routines invoked are:
! \begin{description}
!   \item[LIS\_writeGlobalHeader\_restart](\ref{LIS_writeGlobalHeader_restart})\\
!      writes the global header information
!   \item[LIS\_writeHeader\_restart](\ref{LIS_writeHeader_restart})\\
!      writes the header information for a variable
!   \item[LIS\_closeHeader\_restart](\ref{LIS_closeHeader_restart})\\
!      close the header
!   \item[LIS\_writevar\_restart](\ref{LIS_writevar_restart})\\
!      writes a variable to the restart file
! \end{description}
! 
!EOP 
               
    integer :: l, t 
    real    :: tmptilen(LIS_rc%npatch(n, LIS_rc%lsm_index))
    integer :: dimID(11)
    integer :: sfcrunoff_ID
    integer :: udrrunoff_ID
    integer :: smc_ID
    integer :: sh2o_ID
    integer :: tslb_ID
    integer :: sneqv_ID
    integer :: snowh_ID
    integer :: canwat_ID
    integer :: acsnom_ID
    integer :: acsnow_ID
    integer :: isnow_ID
    integer :: tv_ID
    integer :: tg_ID
    integer :: canice_ID
    integer :: canliq_ID
    integer :: eah_ID
    integer :: tah_ID
    integer :: cm_ID
    integer :: ch_ID
    integer :: fwet_ID
    integer :: sneqvo_ID
    integer :: albold_ID
    integer :: qsnow_ID
    integer :: wslake_ID
    integer :: zwt_ID
    integer :: wa_ID
    integer :: wt_ID
    integer :: tsno_ID
    integer :: zss_ID
    integer :: snowice_ID
    integer :: snowliq_ID
    integer :: lfmass_ID
    integer :: rtmass_ID
    integer :: stmass_ID
    integer :: wood_ID
    integer :: stblcp_ID
    integer :: fastcp_ID
    integer :: lai_ID
    integer :: sai_ID
    integer :: tauss_ID
    integer :: smoiseq_ID
    integer :: smcwtd_ID
    integer :: deeprech_ID
    integer :: rech_ID
    integer :: grain_ID
    integer :: gdd_ID
    integer :: pgs_ID
    integer :: pexp_ID
    integer :: area_ID
    integer :: qrf_ID
    integer :: qspring_ID
    integer :: qslat_ID
    integer :: qrfs_ID
    integer :: qsprings_ID
    integer :: fdepth_ID
    integer :: rivercond_ID
    integer :: riverbed_ID
    integer :: eqzwt_ID
    integer :: irnumsi_ID
    integer :: irnummi_ID
    integer :: irnumfi_ID
    integer :: irwatsi_ID
    integer :: irwatmi_ID
    integer :: irwatfi_ID
    integer :: irsivol_ID
    integer :: irmivol_ID
    integer :: irfivol_ID
    integer :: ireloss_ID
    integer :: irrsplh_ID
    integer :: qtdrain_ID
    integer :: accssoil_ID
    integer :: accqinsur_ID
    integer :: accqseva_ID
    integer :: accetrani_ID
    integer :: accdwater_ID
    integer :: accprcp_ID
    integer :: accecan_ID
    integer :: accetran_ID
    integer :: accedir_ID
    
    ! write the header of the restart file
    call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, &
                                       "NOAHMPnew", &
                                       dim1=NoahMPnew_struc(n)%nsoil+NoahMPnew_struc(n)%nsnow, &
                                       dim2=NoahMPnew_struc(n)%nsoil,                          &
                                       dim3=NoahMPnew_struc(n)%nsnow,                          &
                                       dim4=1,                                                 &
                                       dimID=dimID,                                            &
                                       output_format = trim(wformat))

    ! write the header for state variable sfcrunoff
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sfcrunoff_ID, "SFCRUNOFF", &
                                 "accumulated surface runoff", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable udrrunoff
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, udrrunoff_ID, "UDRRUNOFF", &
                                 "accumulated sub-surface runoff", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable smc
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, smc_ID, "SMC", &
                                 "volumtric soil moisture", &
                                 "m3/m3", vlevels=NoahMPnew_struc(n)%nsoil , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim2") 
   
      ! write the header for state variable sh2o
      !TODO: check dimension of the state variable following "vlevels="
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, sh2o_ID, "SH2O", &
                                   "volumtric liquid soil moisture", &
                                   "m3/m3", vlevels=NoahMPnew_struc(n)%nsoil , valid_min=-99999.0, valid_max=99999.0, &
                                    var_flag = "dim2") 
   
      ! write the header for state variable tslb
      !TODO: check dimension of the state variable following "vlevels="
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, tslb_ID, "TSLB", &
                                   "soil temperature", &
                                   "K", vlevels=NoahMPnew_struc(n)%nsoil , valid_min=-99999.0, valid_max=99999.0, &
                                    var_flag = "dim2")
   
      ! write the header for state variable sneqv
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, sneqv_ID, "SNEQV", &
                                   "snow water equivalent", &
                                   "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable snowh
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, snowh_ID, "SNOWH", &
                                   "physical snow depth", &
                                   "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable canwat
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, canwat_ID, "CANWAT", &
                                   "total canopy water + ice", &
                                   "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable acsnom
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, acsnom_ID, "ACSNOM", &
                                   "accumulated snow melt leaving pack", &
                                   "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable acsnow
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, acsnow_ID, "ACSNOW", &
                                   "accumulated snow on grid", &
                                   "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable isnow
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, isnow_ID, "ISNOW", &
                                   "actual no. of snow layers", &
                                   "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable tv
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, tv_ID, "TV", &
                                   "vegetation leaf temperature", &
                                   "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable tg
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, tg_ID, "TG", &
                                   "bulk ground surface temperature", &
                                   "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable canice
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, canice_ID, "CANICE", &
                                   "canopy-intercepted ice", &
                                   "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable canliq
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, canliq_ID, "CANLIQ", &
                                   "canopy-intercepted liquid water", &
                                   "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable eah
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, eah_ID, "EAH", &
                                   "canopy air vapor pressure", &
                                   "Pa", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable tah
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, tah_ID, "TAH", &
                                   "canopy air temperature", &
                                   "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable cm
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, cm_ID, "CM", &
                                   "bulk momentum drag coefficient", &
                                   "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable ch
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, ch_ID, "CH", &
                                   "bulk sensible heat exchange coefficient", &
                                   "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable fwet
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, fwet_ID, "FWET", &
                                   "wetted or snowed fraction of canopy", &
                                   "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable sneqvo
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, sneqvo_ID, "SNEQVO", &
                                   "snow mass at last time step", &
                                   "mm h2o", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable albold
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, albold_ID, "ALBOLD", &
                                   "snow albedo at last time step", &
                                   "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable qsnow
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, qsnow_ID, "QSNOW", &
                                   "snowfall on the ground", &
                                   "mm/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable wslake
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, wslake_ID, "WSLAKE", &
                                   "lake water storage", &
                                   "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable zwt
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, zwt_ID, "ZWT", &
                                   "water table depth", &
                                   "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable wa
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, wa_ID, "WA", &
                                   "water in aquifer", &
                                   "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable wt
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, wt_ID, "WT", &
                                   "water in aquifer and saturated soil", &
                                   "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable tsno
      !TODO: check dimension of the state variable following "vlevels="
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, tsno_ID, "TSNO", &
                                   "snow layer temperature", &
                                   "K", vlevels=NoahMPnew_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
                                   var_flag = "dim3")
   
      ! write the header for state variable zss
      !TODO: check dimension of the state variable following "vlevels="
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, zss_ID, "ZSS", &
                                   "snow/soil layer depth from snow surface", &
                                   "m", vlevels=NoahMPnew_struc(n)%nsnow+NoahMPnew_struc(n)%nsoil , valid_min=-99999.0, valid_max=99999.0, &
                                   var_flag = "dim1") 
   
      ! write the header for state variable snowice
      !TODO: check dimension of the state variable following "vlevels="
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, snowice_ID, "SNOWICE", &
                                   "snow layer ice", &
                                   "mm", vlevels=NoahMPnew_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
                                    var_flag = "dim3")
   
      ! write the header for state variable snowliq
      !TODO: check dimension of the state variable following "vlevels="
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, snowliq_ID, "SNOWLIQ", &
                                   "snow layer liquid water", &
                                   "mm", vlevels=NoahMPnew_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
                                   var_flag = "dim3")
   
      ! write the header for state variable lfmass
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, lfmass_ID, "LFMASS", &
                                   "leaf mass", &
                                   "g/m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable rtmass
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, rtmass_ID, "RTMASS", &
                                   "mass of fine roots", &
                                   "g/m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable stmass
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, stmass_ID, "STMASS", &
                                   "stem mass", &
                                   "g/m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable wood
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, wood_ID, "WOOD", &
                                   "mass of wood (including woody roots)", &
                                   "g/m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable stblcp
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, stblcp_ID, "STBLCP", &
                                   "stable carbon in deep soil", &
                                   "g/m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
      ! write the header for state variable fastcp
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, fastcp_ID, "FASTCP", &
                                   "short-lived carbon in shallow soil", &
                                   "g/m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
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
      ! write the header for state variable tauss
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, tauss_ID, "TAUSS", &
                                   "snow age factor", &
                                   "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

      ! for MMF groundwater
      if (NoahMPnew_struc(n)%runsub_opt == 5) then 
      ! write the header for state variable smoiseq
      !TODO: check dimension of the state variable following "vlevels="
      !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
      call LIS_writeHeader_restart(ftn, n, dimID, smoiseq_ID, "SMOISEQ", &
                                   "equilibrium volumetric soil moisture content", &
                                   "m3/m3", vlevels=NoahMPnew_struc(n)%nsoil , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim2") 
 
    ! write the header for state variable smcwtd
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, smcwtd_ID, "SMCWTD", &
                                 "soil moisture content in the layer to the water table when deep", &
                                 "m3/m3", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable deeprech
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, deeprech_ID, "DEEPRECH", &
                                 "recharge to the water table when deep", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable rech
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, rech_ID, "RECH", &
                                 "recharge to the water table (diagnostic)", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable pexp
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, pexp_ID, "PEXP", &
                                 "groundwater expotential parameter", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable area
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, area_ID, "AREA", &
                                 "river area", &
                                 "m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable qrf
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qrf_ID, "QRF", &
                                 "groundwater baselow", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable qspring
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qspring_ID, "QSPRING", &
                                 "seeping water", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable qslat
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qslat_ID, "QSLAT", &
                                 "accumulated lateral flow", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable qrfs
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qrfs_ID, "QRFS", &
                                 "accumulated GW baseflow", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable qsprings
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qsprings_ID, "QSPRINGS", &
                                 "accumulated seeping water", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable fdepth
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, fdepth_ID, "FDEPTH", &
                                 "depth", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable rivercond
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, rivercond_ID, "RIVERCOND", &
                                 "river conductivity", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable riverbed
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, riverbed_ID, "RIVERBED", &
                                 "riverbed depth", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable eqzwt
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, eqzwt_ID, "EQZWT", &
                                 "equilibrium water table depth", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    endif

    ! for irrigation
    if (NoahMPnew_struc(n)%irr_opt >0) then
    ! write the header for state variable irnumsi
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, irnumsi_ID, "IRNUMSI", &
                                 "sprinkler irrigation count", &
                                 "-", vlevels=1, valid_min=-99999, valid_max=99999)
    ! write the header for state variable irnummi
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, irnummi_ID, "IRNUMMI", &
                                 "micro irrigation count", &
                                 "-", vlevels=1, valid_min=-99999, valid_max=99999)
    ! write the header for state variable irnumfi
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, irnumfi_ID, "IRNUMFI", &
                                 "flood irrigation count", &
                                 "-", vlevels=1, valid_min=-99999, valid_max=99999)
    ! write the header for state variable irwatsi
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, irwatsi_ID, "IRWATSI", &
                                 "sprinkler irrigation water amount", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable irwatmi
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, irwatmi_ID, "IRWATMI", &
                                 "micro irrigation water amount", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable irwatfi
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, irwatfi_ID, "IRWATFI", &
                                 "flood irrigation water amount", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable irsivol
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, irsivol_ID, "IRSIVOL", &
                                 "sprinkler irrigation water volume", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable irmivol
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, irmivol_ID, "IRMIVOL", &
                                 "micro irrigation water volume", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable irfivol
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, irfivol_ID, "IRFIVOL", &
                                 "flood irrigation water volume", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable ireloss
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, ireloss_ID, "IRELOSS", &
                                 "loss of irrigation water to evaporation", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable irrsplh
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, irrsplh_ID, "IRRSPLH", &
                                 "latent heating from sprinkler evaporation", &
                                 "W/m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    endif

    ! for tile drainage
    if (NoahMPnew_struc(n)%tdrn_opt >0) then
    ! write the header for state variable qtdrain
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qtdrain_ID, "QTDRAIN", &
                                 "accumulated tile drainage discharge", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    endif

    ! write the header for state variable grain
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, grain_ID, "GRAIN", &
                                 "mass of grain XING", &
                                 "g/m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable gdd
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, gdd_ID, "GDD", &
                                 "growing degree days XING (based on 10C)", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable pgs
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, pgs_ID, "PGS", &
                                 "growing degree days XING", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! for additional restart variables
    ! write the header for state variable accssoil
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, accssoil_ID, "ACC_SSOIL", &
                                 "accumulated ground heat flux", &
                                 "W/m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable accqinsur
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, accqinsur_ID, "ACC_QINSUR", &
                                 "accumulated soil surface water flux", &
                                 "m/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable accqseva
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, accqseva_ID, "ACC_QSEVA", &
                                 "accumulated soil surface evaporation", &
                                 "m/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable accetrani
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, accetrani_ID, "ACC_ETRANI", &
                                 "accumulated plant transpiration", &
                                 "m/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable accdwater
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, accdwater_ID, "ACC_DWATER", &
                                 "accumulated water storage change", &
                                 "m/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable accprcp
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, accprcp_ID, "ACC_PRCP", &
                                 "accumulated precipitation", &
                                 "m/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable accecan
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, accecan_ID, "ACC_ECAN", &
                                 "accumulated canopy evaporation", &
                                 "m/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable accetran
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, accetran_ID, "ACC_ETRAN", &
                                 "accumulated transpiration", &
                                 "m/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable accedir
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, accedir_ID, "ACC_EDIR", &
                                 "accumulated net soil evaporation", &
                                 "m/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

 
    ! close header of restart file
    call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, NoahMPnew_struc(n)%rstInterval)

    ! write state variables into restart file
    ! accumulated surface runoff
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%sfcrunoff, &
                              varid=sfcrunoff_ID, dim=1, wformat=wformat)

    ! accumulated sub-surface runoff
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%udrrunoff, &
                              varid=udrrunoff_ID, dim=1, wformat=wformat)

    ! volumtric soil moisture
    do l=1, NoahMPnew_struc(n)%nsoil   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = NoahMPnew_struc(n)%noahmpnew(t)%smc(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=smc_ID, dim=l, wformat=wformat)
    enddo
    ! volumtric liquid soil moisture
    do l=1, NoahMPnew_struc(n)%nsoil   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = NoahMPnew_struc(n)%noahmpnew(t)%sh2o(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=sh2o_ID, dim=l, wformat=wformat)
    enddo
    ! soil temperature
    do l=1, NoahMPnew_struc(n)%nsoil   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = NoahMPnew_struc(n)%noahmpnew(t)%tslb(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=tslb_ID, dim=l, wformat=wformat)
    enddo
    ! snow water equivalent
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%sneqv, &
                              varid=sneqv_ID, dim=1, wformat=wformat)

    ! physical snow depth
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%snowh, &
                              varid=snowh_ID, dim=1, wformat=wformat)

    ! total canopy water + ice
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%canwat, &
                              varid=canwat_ID, dim=1, wformat=wformat)

    ! accumulated snow melt leaving pack
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%acsnom, &
                              varid=acsnom_ID, dim=1, wformat=wformat)

    ! accumulated snow on grid
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%acsnow, &
                              varid=acsnow_ID, dim=1, wformat=wformat)

    ! actual no. of snow layers
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%isnow, &
                              varid=isnow_ID, dim=1, wformat=wformat)

    ! vegetation leaf temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%tv, &
                              varid=tv_ID, dim=1, wformat=wformat)

    ! bulk ground surface temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%tg, &
                              varid=tg_ID, dim=1, wformat=wformat)

    ! canopy-intercepted ice
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%canice, &
                              varid=canice_ID, dim=1, wformat=wformat)

    ! canopy-intercepted liquid water
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%canliq, &
                              varid=canliq_ID, dim=1, wformat=wformat)

    ! canopy air vapor pressure
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%eah, &
                              varid=eah_ID, dim=1, wformat=wformat)

    ! canopy air temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%tah, &
                              varid=tah_ID, dim=1, wformat=wformat)

    ! bulk momentum drag coefficient
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%cm, &
                              varid=cm_ID, dim=1, wformat=wformat)

    ! bulk sensible heat exchange coefficient
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%ch, &
                              varid=ch_ID, dim=1, wformat=wformat)

    ! wetted or snowed fraction of canopy
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%fwet, &
                              varid=fwet_ID, dim=1, wformat=wformat)

    ! snow mass at last time step
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%sneqvo, &
                              varid=sneqvo_ID, dim=1, wformat=wformat)

    ! snow albedo at last time step
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%albold, &
                              varid=albold_ID, dim=1, wformat=wformat)

    ! snowfall on the ground
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%qsnow, &
                              varid=qsnow_ID, dim=1, wformat=wformat)

    ! lake water storage
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%wslake, &
                              varid=wslake_ID, dim=1, wformat=wformat)

    ! water table depth
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%zwt, &
                              varid=zwt_ID, dim=1, wformat=wformat)

    ! water in the "aquifer"
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%wa, &
                              varid=wa_ID, dim=1, wformat=wformat)

    ! water in aquifer and saturated soil
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%wt, &
                              varid=wt_ID, dim=1, wformat=wformat)

    ! snow layer temperature
    do l=1, NoahMPnew_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = NoahMPnew_struc(n)%noahmpnew(t)%tsno(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=tsno_ID, dim=l, wformat=wformat)
    enddo
    ! snow/soil layer depth from snow surface
    do l=1, NoahMPnew_struc(n)%nsnow+NoahMPnew_struc(n)%nsoil   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = NoahMPnew_struc(n)%noahmpnew(t)%zss(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=zss_ID, dim=l, wformat=wformat)
    enddo
    ! snow layer ice
    do l=1, NoahMPnew_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = NoahMPnew_struc(n)%noahmpnew(t)%snowice(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=snowice_ID, dim=l, wformat=wformat)
    enddo
    ! snow layer liquid water
    do l=1, NoahMPnew_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = NoahMPnew_struc(n)%noahmpnew(t)%snowliq(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=snowliq_ID, dim=l, wformat=wformat)
    enddo
    ! leaf mass
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%lfmass, &
                              varid=lfmass_ID, dim=1, wformat=wformat)

    ! mass of fine roots
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%rtmass, &
                              varid=rtmass_ID, dim=1, wformat=wformat)

    ! stem mass
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%stmass, &
                              varid=stmass_ID, dim=1, wformat=wformat)

    ! mass of wood (including woody roots)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%wood, &
                              varid=wood_ID, dim=1, wformat=wformat)

    ! stable carbon in deep soil
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%stblcp, &
                              varid=stblcp_ID, dim=1, wformat=wformat)

    ! short-lived carbon in shallow soil
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%fastcp, &
                              varid=fastcp_ID, dim=1, wformat=wformat)

    ! leaf area index
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%lai, &
                              varid=lai_ID, dim=1, wformat=wformat)

    ! stem area index
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%sai, &
                              varid=sai_ID, dim=1, wformat=wformat)

    ! snow age factor
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%tauss, &
                              varid=tauss_ID, dim=1, wformat=wformat)

    ! equilibrium volumetric soil moisture content
    do l=1, NoahMPnew_struc(n)%nsoil   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = NoahMPnew_struc(n)%noahmpnew(t)%smoiseq(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=smoiseq_ID, dim=l, wformat=wformat)
    enddo
    ! soil moisture content in the layer to the water table when deep
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%smcwtd, &
                              varid=smcwtd_ID, dim=1, wformat=wformat)

    ! recharge to the water table when deep
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%deeprech, &
                              varid=deeprech_ID, dim=1, wformat=wformat)

    ! recharge to the water table (diagnostic)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%rech, &
                              varid=rech_ID, dim=1, wformat=wformat)

    ! mass of grain XING
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%grain, &
                              varid=grain_ID, dim=1, wformat=wformat)

    ! growing degree days XING (based on 10C)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%gdd, &
                              varid=gdd_ID, dim=1, wformat=wformat)

    ! growing degree days XING
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, NoahMPnew_struc(n)%noahmpnew%pgs, &
                              varid=pgs_ID, dim=1, wformat=wformat)

    ! optional gecros crop
!    do l=1, 60  ! TODO: check loop
!        tmptilen = 0
!        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
!            tmptilen(t) = NoahMPnew_struc(n)%noahmpnew(t)%gecros_state(l)
!        enddo
!        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
!                                  varid=gecros_state_ID, dim=l, wformat=wformat)
!    enddo
end subroutine NoahMPnew_dump_restart
