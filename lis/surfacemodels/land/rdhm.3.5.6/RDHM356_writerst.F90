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
! !ROUTINE: RDHM356_writerst
! \label{RDHM356_writerst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   11/5/13: Shugong Wang; initial implementation for LIS 7 and RDHM356
!
! !INTERFACE:
subroutine RDHM356_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use RDHM356_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for RDHM356.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[RDHM356\_dump\_restart](\ref{RDHM356_dump_restart}) \newline
!   writes the RDHM356 variables into the restart file
! \end{description}
!EOP

    character(len=LIS_CONST_PATH_LEN) :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn
    integer       :: status
    
    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "RDHM356 restart alarm")
    
    ! set restart file format (read from LIS configration file_
    wformat = trim(RDHM356_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", "RDHM356",&
                                             wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status,"Error in nf90_open in RDHM356_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, "Error in nf90_open in RDHM356_writerst")
#endif
             endif
        endif
    
        call RDHM356_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in RDHM356_writerst")
#endif
            endif
            write(LIS_logunit, *) "RDHM356 archive restart written: ", trim(filen)
        endif
    endif
end subroutine RDHM356_writerst

!BOP
!
! !ROUTINE: RDHM356_dump_restart
! \label{RDHM356_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  11/5/13: Shugong Wang, initial implementation for LIS 7 and RDHM356
! !INTERFACE:
subroutine RDHM356_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use RDHM356_lsmMod

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
!  The following is the list of variables written in the RDHM356
!  restart file:
!\small
!  \begin{verbatim}
!    nc, nr, ntiles  - grid and tile space dimensions
!    UZTWC           - RDHM356 upper zone tension water storage content [mm]
!    UZFWC           - RDHM356 upper zone free water storage content [mm]
!    LZTWC           - RDHM356 lower zone tension water storage content [mm]
!    LZFPC           - RDHM356 lower zone primary free water storage content [mm]
!    LZFSC           - RDHM356 lower zone supplemental free water storage content [mm]
!    ADIMC           - RDHM356 additional impervious area content [mm]
!    TS0             - RDHM356 first soil layer temperature [C]
!    TS1             - RDHM356 second soil layer temperature [C]
!    TS2             - RDHM356 third soil layer temperature [C]
!    TS3             - RDHM356 fourth soil layer temperature [C]
!    TS4             - RDHM356 fifth soil layer temperature [C]
!    UZTWH           - RDHM356 unfrozen upper zone tension water [mm]
!    UZFWH           - RDHM356 unfrozen uppeer zone free water [mm]
!    LZTWH           - RDHM356 unfrozen lower zone tension water [mm]
!    LZFSH           - RDHM356 unfrozen lower zone supplemental free water [mm]
!    LZFPH           - RDHM356 unfrozen lower zone primary free water [mm]
!    SMC             - RDHM356 volumetric content of total soil moisture at each layer [m^3 m-3]
!    SH2O            - RDHM356 volumetric content of liquid soil moisture at each layer [m^3 m-3]
!    WE              - RDHM356 snow water equivalent without liquid water [mm]
!    LIQW            - RDHM356 liquid water in snow [mm]
!    NEGHS           - RDHM356 negative snow heat [mm]
!    TINDEX          - RDHM356 antecedent temperature index [C]
!    ACCMAX          - RDHM356 cumulated snow water including liquid [mm]
!    SNDPT           - RDHM356 snow depth [cm]
!    SNTMP           - RDHM356 average snow temperature [C]
!    SB              - RDHM356 the last highest snow water equivalent before any snow fall [C]
!    SBAESC          - RDHM356 internal snow state during melt & new snow fall (checked with Victor) [-]
!    SBWS            - RDHM356 internal snow state during melt & new snow fall (checked with Victor) [-]
!    STORAGE         - RDHM356 snow liquid water attenuation storage [mm]
!    AEADJ           - RDHM356 adjusted areal snow cover fraction [-]
!    EXLAG           - RDHM356 array of lagged liquid water values [-]
!    NEXLAG          - RDHM356 number of ordinates in lagged liquid water array (EXLAG) [-]
!    TA_PREV         - RDHM356 air temperature of previous time step [-]
!  \end{verbatim}
!\normalsize
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
    integer :: UZTWC_ID
    integer :: UZFWC_ID
    integer :: LZTWC_ID
    integer :: LZFPC_ID
    integer :: LZFSC_ID
    integer :: ADIMC_ID
    integer :: TS0_ID
    integer :: TS1_ID
    integer :: TS2_ID
    integer :: TS3_ID
    integer :: TS4_ID
    integer :: UZTWH_ID
    integer :: UZFWH_ID
    integer :: LZTWH_ID
    integer :: LZFSH_ID
    integer :: LZFPH_ID
    integer :: SMC_ID
    integer :: SH2O_ID
    integer :: WE_ID
    integer :: LIQW_ID
    integer :: NEGHS_ID
    integer :: TINDEX_ID
    integer :: ACCMAX_ID
    integer :: SNDPT_ID
    integer :: SNTMP_ID
    integer :: SB_ID
    integer :: SBAESC_ID
    integer :: SBWS_ID
    integer :: STORAGE_ID
    integer :: AEADJ_ID
    integer :: EXLAG_ID
    integer :: NEXLAG_ID
    integer :: TA_PREV_ID
    integer :: CH_ID, CM_ID
    
    ! write the header of the restart file
    call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, &
                                       "RDHM356", dim1=6, dim2=7, dimID=dimID, &
                                       output_format = trim(wformat))

    ! write the header for state variable UZTWC
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, UZTWC_ID, "UZTWC", &
                                 "upper zone tension water storage content", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable UZFWC
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, UZFWC_ID, "UZFWC", &
                                 "upper zone free water storage content", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable LZTWC
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, LZTWC_ID, "LZTWC", &
                                 "lower zone tension water storage content", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable LZFPC
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, LZFPC_ID, "LZFPC", &
                                 "lower zone primary free water storage content", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable LZFSC
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, LZFSC_ID, "LZFSC", &
                                 "lower zone supplemental free water storage content", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable ADIMC
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, ADIMC_ID, "ADIMC", &
                                 "additional impervious area content", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable TS0
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, TS0_ID, "TS0", &
                                 "first soil layer temperature", &
                                 "C", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable TS1
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, TS1_ID, "TS1", &
                                 "second soil layer temperature", &
                                 "C", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable TS2
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, TS2_ID, "TS2", &
                                 "third soil layer temperature", &
                                 "C", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable TS3
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, TS3_ID, "TS3", &
                                 "fourth soil layer temperature", &
                                 "C", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable TS4
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, TS4_ID, "TS4", &
                                 "fifth soil layer temperature", &
                                 "C", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable UZTWH
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, UZTWH_ID, "UZTWH", &
                                 "unfrozen upper zone tension water", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable UZFWH
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, UZFWH_ID, "UZFWH", &
                                 "unfrozen uppeer zone free water", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable LZTWH
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, LZTWH_ID, "LZTWH", &
                                 "unfrozen lower zone tension water", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable LZFSH
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, LZFSH_ID, "LZFSH", &
                                 "unfrozen lower zone supplemental free water", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable LZFPH
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, LZFPH_ID, "LZFPH", &
                                 "unfrozen lower zone primary free water", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable SMC
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SMC_ID, "SMC", &
                                 "volumetric content of total soil moisture at each layer", &
                                 "m^3 m-3", vlevels=6, valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")
 
    ! write the header for state variable SH2O
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SH2O_ID, "SH2O", &
                                 "volumetric content of liquid soil moisture at each layer", &
                                 "m^3 m-3", vlevels=6, valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")
 
    ! write the header for state variable WE
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, WE_ID, "WE", &
                                 "snow water equivalent without liquid water", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable LIQW
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, LIQW_ID, "LIQW", &
                                 "liquid water in snow", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable NEGHS
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, NEGHS_ID, "NEGHS", &
                                 "negative snow heat", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable TINDEX
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, TINDEX_ID, "TINDEX", &
                                 "antecedent temperature index", &
                                 "C", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable ACCMAX
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, ACCMAX_ID, "ACCMAX", &
                                 "cumulated snow water including liquid", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable SNDPT
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNDPT_ID, "SNDPT", &
                                 "snow depth", &
                                 "cm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable SNTMP
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNTMP_ID, "SNTMP", &
                                 "average snow temperature", &
                                 "C", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable SB
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SB_ID, "SB", &
                                 "the last highest snow water equivalent before any snow fall", &
                                 "C", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable SBAESC
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SBAESC_ID, "SBAESC", &
                                 "internal snow state during melt & new snow fall (checked with Victor)", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable SBWS
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SBWS_ID, "SBWS", &
                                 "internal snow state during melt & new snow fall (checked with Victor)", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable STORAGE
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, STORAGE_ID, "STORAGE", &
                                 "snow liquid water attenuation storage", &
                                 "mm", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable AEADJ
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, AEADJ_ID, "AEADJ", &
                                 "adjusted areal snow cover fraction", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! write the header for state variable EXLAG
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, EXLAG_ID, "EXLAG", &
                                 "array of lagged liquid water values", &
                                 "-", vlevels=7, valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim2")
 
    ! write the header for state variable NEXLAG
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, NEXLAG_ID, "NEXLAG", &
                                 "number of ordinates in lagged liquid water array (EXLAG)", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! CH
    call LIS_writeHeader_restart(ftn, n, dimID, CH_ID, "CH", &
                                 "coefficient of heat exchange", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! CM
    call LIS_writeHeader_restart(ftn, n, dimID, CM_ID, "CM", &
                                 "coefficient of momentum exchange", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    
    ! write the header for state variable TA_PREV
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, TA_PREV_ID, "TA_PREV", &
                                 "air temperature of previous time step", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
 
    ! close header of restart file
    call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, RDHM356_struc(n)%rstInterval)

    ! write state variables into restart file
    ! upper zone tension water storage content
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%UZTWC, &
                              varid=UZTWC_ID, dim=1, wformat=wformat)

    ! upper zone free water storage content
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%UZFWC, &
                              varid=UZFWC_ID, dim=1, wformat=wformat)

    ! lower zone tension water storage content
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZTWC, &
                              varid=LZTWC_ID, dim=1, wformat=wformat)

    ! lower zone primary free water storage content
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZFPC, &
                              varid=LZFPC_ID, dim=1, wformat=wformat)

    ! lower zone supplemental free water storage content
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZFSC, &
                              varid=LZFSC_ID, dim=1, wformat=wformat)

    ! additional impervious area content
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%ADIMC, &
                              varid=ADIMC_ID, dim=1, wformat=wformat)

    ! first soil layer temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TS0, &
                              varid=TS0_ID, dim=1, wformat=wformat)

    ! second soil layer temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TS1, &
                              varid=TS1_ID, dim=1, wformat=wformat)

    ! third soil layer temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TS2, &
                              varid=TS2_ID, dim=1, wformat=wformat)

    ! fourth soil layer temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TS3, &
                              varid=TS3_ID, dim=1, wformat=wformat)

    ! fifth soil layer temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TS4, &
                              varid=TS4_ID, dim=1, wformat=wformat)

    ! unfrozen upper zone tension water
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%UZTWH, &
                              varid=UZTWH_ID, dim=1, wformat=wformat)

    ! unfrozen uppeer zone free water
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%UZFWH, &
                              varid=UZFWH_ID, dim=1, wformat=wformat)

    ! unfrozen lower zone tension water
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZTWH, &
                              varid=LZTWH_ID, dim=1, wformat=wformat)

    ! unfrozen lower zone supplemental free water
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZFSH, &
                              varid=LZFSH_ID, dim=1, wformat=wformat)

    ! unfrozen lower zone primary free water
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LZFPH, &
                              varid=LZFPH_ID, dim=1, wformat=wformat)

    ! volumetric content of total soil moisture at each layer
    do l=1, 6  ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = RDHM356_struc(n)%rdhm356(t)%SMC(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SMC_ID, dim=l, wformat=wformat)
    enddo
    ! volumetric content of liquid soil moisture at each layer
    do l=1, 6  ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = RDHM356_struc(n)%rdhm356(t)%SH2O(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SH2O_ID, dim=l, wformat=wformat)
    enddo
    ! snow water equivalent without liquid water
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%WE, &
                              varid=WE_ID, dim=1, wformat=wformat)

    ! liquid water in snow
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%LIQW, &
                              varid=LIQW_ID, dim=1, wformat=wformat)

    ! negative snow heat
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%NEGHS, &
                              varid=NEGHS_ID, dim=1, wformat=wformat)

    ! antecedent temperature index
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TINDEX, &
                              varid=TINDEX_ID, dim=1, wformat=wformat)

    ! cumulated snow water including liquid
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%ACCMAX, &
                              varid=ACCMAX_ID, dim=1, wformat=wformat)

    ! snow depth
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%SNDPT, &
                              varid=SNDPT_ID, dim=1, wformat=wformat)

    ! average snow temperature
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%SNTMP, &
                              varid=SNTMP_ID, dim=1, wformat=wformat)

    ! the last highest snow water equivalent before any snow fall
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%SB, &
                              varid=SB_ID, dim=1, wformat=wformat)

    ! internal snow state during melt & new snow fall (checked with Victor)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%SBAESC, &
                              varid=SBAESC_ID, dim=1, wformat=wformat)

    ! internal snow state during melt & new snow fall (checked with Victor)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%SBWS, &
                              varid=SBWS_ID, dim=1, wformat=wformat)

    ! snow liquid water attenuation storage
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%STORAGE, &
                              varid=STORAGE_ID, dim=1, wformat=wformat)

    ! adjusted areal snow cover fraction
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%AEADJ, &
                              varid=AEADJ_ID, dim=1, wformat=wformat)

    ! array of lagged liquid water values
    do l=1, 7  ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = RDHM356_struc(n)%rdhm356(t)%EXLAG(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=EXLAG_ID, dim=l, wformat=wformat)
    enddo
    ! number of ordinates in lagged liquid water array (EXLAG)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%NEXLAG, &
                              varid=NEXLAG_ID, dim=1, wformat=wformat)

    ! air temperature of previous time step
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%TA_PREV, &
                              varid=TA_PREV_ID, dim=1, wformat=wformat)
    ! ch
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%CH, &
                              varid=CH_ID, dim=1, wformat=wformat)
    ! cm
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RDHM356_struc(n)%rdhm356%CM, &
                              varid=CM_ID, dim=1, wformat=wformat)

end subroutine RDHM356_dump_restart
