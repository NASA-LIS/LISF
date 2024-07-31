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
! !ROUTINE: Crocus81_writerst
! \label{Crocus81_writerst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/18/19: Mahdi Navari, Shugong Wang; initial implementation for LIS 7 and Crocus81
!
! !INTERFACE:
subroutine Crocus81_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use Crocus81_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for Crocus81.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory})\\
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename})\\
!  generates a timestamped restart filename
! \item[Crocus81\_dump\_restart](\ref{Crocus81_dump_restart})\\
!   writes the Crocus81 variables into the restart file
! \end{description}
!EOP

    character(len=LIS_CONST_PATH_LEN) :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn
    integer       :: status
    character*3        :: fnest  ! MN added 

    ! set restart alarm
    ! alarmCheck = LIS_isAlarmRinging(LIS_rc, "Crocus81 restart alarm")
    write(fnest,'(i3.3)') n   ! MN added   Bug in the toolkit
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "CROCUS81 restart alarm "//trim(fnest)) ! MN added  Bug in the toolkit
    ! set restart file format (read from LIS configration file_
    wformat = trim(CROCUS81_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", "CROCUS81",&
                                             wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status,"Error in nf90_open in Crocus81_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, "Error in nf90_open in Crocus81_writerst")
#endif
             endif
        endif
    
        call Crocus81_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in Crocus81_writerst")
#endif
            endif
            write(LIS_logunit, *) "Crocus81 archive restart written: ", trim(filen)
        endif
    endif
end subroutine Crocus81_writerst

!BOP
!
! !ROUTINE: Crocus81_dump_restart
! \label{Crocus81_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  10/18/19: Mahdi Navari, Shugong Wang, initial implementation for LIS 7 and Crocus81
! !INTERFACE:
subroutine Crocus81_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use Crocus81_lsmMod

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
!  The following is the list of variables written in the Crocus81
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
    integer :: dimID(11)  !6
    integer :: SNOWSWE_ID
    integer :: SNOWRHO_ID
    integer :: SNOWHEAT_ID
    integer :: SNOWALB_ID
    integer :: SNOWGRAN1_ID
    integer :: SNOWGRAN2_ID
    integer :: SNOWHIST_ID
    integer :: SNOWAGE_ID
    integer :: SNOWLIQ_ID
    integer :: SNOWTEMP_ID
    integer :: SNOWDZ_ID
    integer :: GRNDFLUX_ID
    integer :: SNDRIFT_ID
    integer :: RI_n_ID
    integer :: CDSNOW_ID
    integer :: USTARSNOW_ID
    integer :: CHSNOW_ID
    integer :: SNOWMAK_dz_ID
    
    ! write the header of the restart file
   call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, &
                                       "CROCUS", &
                                       dim1=CROCUS81_struc(n)%nsnow, &
                                       dim2=CROCUS81_struc(n)%nimpur, &
                                       dim3=3,   &  ! reserved for spectral albedo for case B92 for case T17 the dim will be 186
                                       dimID=dimID, &
                                       output_format = trim(wformat))

    ! write the header for state variable SNOWSWE
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWSWE_ID, "SNOWSWE", &
                                 "Snow layer(s) liquid Water Equivalent (SWE:kg m-2)", &
                                 "kg/m2", vlevels=CROCUS81_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable SNOWRHO
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWRHO_ID, "SNOWRHO", &
                                 "Snow layer(s) averaged density (kg/m3)", &
                                 "kg/m3", vlevels=CROCUS81_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable SNOWHEAT
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWHEAT_ID, "SNOWHEAT", &
                                 "Snow layer(s) heat content (J/m2)", &
                                 "J/m2", vlevels=CROCUS81_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable SNOWALB
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWALB_ID, "SNOWALB", &
                                 "snow surface albedo", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SNOWGRAN1
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWGRAN1_ID, "SNOWGRAN1", &
                                 "Snow layers grain feature 1", &
                                 "-", vlevels=CROCUS81_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable SNOWGRAN2
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWGRAN2_ID, "SNOWGRAN2", &
                                 "Snow layer grain feature 2", &
                                 "-", vlevels=CROCUS81_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable SNOWHIST
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWHIST_ID, "SNOWHIST", &
                                 "Snow layer grain historical parameter (only for non dendritic snow) (-) in {0-5}", &
                                 "-", vlevels=CROCUS81_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable SNOWAGE
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWAGE_ID, "SNOWAGE", &
                                 "Age since snowfall (day)", &
                                 "day", vlevels=CROCUS81_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable SNOWLIQ
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWLIQ_ID, "SNOWLIQ", &
                                 "Snow layer(s) liquid water content (m)", &
                                 "m", vlevels=CROCUS81_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable SNOWTEMP
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWTEMP_ID, "SNOWTEMP", &
                                 "Snow layer(s) temperature (K)", &
                                 "K", vlevels=CROCUS81_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable SNOWDZ
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWDZ_ID, "SNOWDZ", &
                                 "Snow layer(s) thickness (m)", &
                                 "m", vlevels=CROCUS81_struc(n)%nsnow , valid_min=-99999.0, valid_max=99999.0, &
var_flag = "dim1") !TODO: replace "xxxx" with correct "dimx"
 
    ! write the header for state variable GRNDFLUX
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, GRNDFLUX_ID, "GRNDFLUX", &
                                 "Soil/snow interface heat flux (W/m2)", &
                                 "W/m2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SNDRIFT
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNDRIFT_ID, "SNDRIFT", &
                                 "Blowing snow sublimation (kg/m2/s) NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592)", &
                                 "kg/m2/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable RI_n
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, RI_n_ID, "RI_N", &
                                 "Richardson number (-)  NOTE: RI has not been initialized in CALL_MODEL (If not OMED initalized to undefined in the snow3L_isba.F90)", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable CDSNOW
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, CDSNOW_ID, "CDSNOW", &
                                 "Drag coefficient for momentum over snow (-)", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable USTARSNOW
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, USTARSNOW_ID, "USTARSNOW", &
                                 "Friction velocity over snow (m/s);", &
                                 "m/s", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable CHSNOW
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, CHSNOW_ID, "CHSNOW", &
                                 "Drag coefficient for heat over snow  (-)", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable SNOWMAK_dz
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, SNOWMAK_dz_ID, "SNOWMAK_DZ", &
                                 "Snowmaking thickness (m)", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! close header of restart file
    call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, CROCUS81_struc(n)%rstInterval)

    ! write state variables into restart file
    ! Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
    do l=1, CROCUS81_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = CROCUS81_struc(n)%crocus81(t)%SNOWSWE(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SNOWSWE_ID, dim=l, wformat=wformat)
    enddo
    ! Snow layer(s) averaged density (kg/m3)
    do l=1, CROCUS81_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = CROCUS81_struc(n)%crocus81(t)%SNOWRHO(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SNOWRHO_ID, dim=l, wformat=wformat)
    enddo
    ! Snow layer(s) heat content (J/m2)
    do l=1, CROCUS81_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = CROCUS81_struc(n)%crocus81(t)%SNOWHEAT(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SNOWHEAT_ID, dim=l, wformat=wformat)
    enddo
    ! snow surface albedo
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%SNOWALB, &
                              varid=SNOWALB_ID, dim=1, wformat=wformat)

    ! Snow layers grain feature 1
    do l=1, CROCUS81_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = CROCUS81_struc(n)%crocus81(t)%SNOWGRAN1(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SNOWGRAN1_ID, dim=l, wformat=wformat)
    enddo
    ! Snow layer grain feature 2
    do l=1, CROCUS81_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SNOWGRAN2_ID, dim=l, wformat=wformat)
    enddo
    ! Snow layer grain historical parameter (only for non dendritic snow) (-) in {0-5}
    do l=1, CROCUS81_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = CROCUS81_struc(n)%crocus81(t)%SNOWHIST(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SNOWHIST_ID, dim=l, wformat=wformat)
    enddo
    ! Age since snowfall (day)
    do l=1, CROCUS81_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = CROCUS81_struc(n)%crocus81(t)%SNOWAGE(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SNOWAGE_ID, dim=l, wformat=wformat)
    enddo
    ! Snow layer(s) liquid water content (m)
    do l=1, CROCUS81_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = CROCUS81_struc(n)%crocus81(t)%SNOWLIQ(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SNOWLIQ_ID, dim=l, wformat=wformat)
    enddo
    ! Snow layer(s) temperature (K)
    do l=1, CROCUS81_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = CROCUS81_struc(n)%crocus81(t)%SNOWTEMP(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SNOWTEMP_ID, dim=l, wformat=wformat)
    enddo
    ! Snow layer(s) thickness (m)
    do l=1, CROCUS81_struc(n)%nsnow   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = CROCUS81_struc(n)%crocus81(t)%SNOWDZ(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=SNOWDZ_ID, dim=l, wformat=wformat)
    enddo
    ! Soil/snow interface heat flux (W/m2)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%GRNDFLUX, &
                              varid=GRNDFLUX_ID, dim=1, wformat=wformat)

    ! Blowing snow sublimation (kg/m2/s) NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%SNDRIFT, &
                              varid=SNDRIFT_ID, dim=1, wformat=wformat)

    ! Richardson number (-)  NOTE: RI has not been initialized in CALL_MODEL (If not OMED initalized to undefined in the snow3L_isba.F90)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%RI_n, &
                              varid=RI_n_ID, dim=1, wformat=wformat)

    ! Drag coefficient for momentum over snow (-)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%CDSNOW, &
                              varid=CDSNOW_ID, dim=1, wformat=wformat)

    ! Friction velocity over snow (m/s);
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%USTARSNOW, &
                              varid=USTARSNOW_ID, dim=1, wformat=wformat)

    ! Drag coefficient for heat over snow  (-)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%CHSNOW, &
                              varid=CHSNOW_ID, dim=1, wformat=wformat)

    ! Snowmaking thickness (m)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, CROCUS81_struc(n)%crocus81%SNOWMAK_dz, &
                              varid=SNOWMAK_dz_ID, dim=1, wformat=wformat)

end subroutine Crocus81_dump_restart
