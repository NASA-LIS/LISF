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
! !ROUTINE: AWRAL600_writerst
! \label{AWRAL600_writerst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   12/18/18: Wendy Sharples, Shugong Wang; initial implementation for LIS 7 and AWRAL600
!
! !INTERFACE:
subroutine AWRAL600_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use AWRAL600_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for AWRAL600.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory})\\
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename})\\
!  generates a timestamped restart filename
! \item[AWRAL600\_dump\_restart](\ref{AWRAL600_dump_restart})\\
!   writes the AWRAL600 variables into the restart file
! \end{description}
!EOP

    character(len=LIS_CONST_PATH_LEN) :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn
    integer       :: status
    
    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "AWRAL600 restart alarm")
    
    ! set restart file format (read from LIS configration file_
    wformat = trim(AWRAL600_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", "AWRAL600",&
                                             wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status,"Error in nf90_open in AWRAL600_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, "Error in nf90_open in AWRAL600_writerst")
#endif
             endif
        endif
    
        call AWRAL600_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in AWRAL600_writerst")
#endif
            endif
            write(LIS_logunit, *) "[INFO] AWRAL600 archive restart written: ", trim(filen)
        endif
    endif
end subroutine AWRAL600_writerst

!BOP
!
! !ROUTINE: AWRAL600_dump_restart
! \label{AWRAL600_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  12/18/18: Wendy Sharples, Shugong Wang, initial implementation for LIS 7 and AWRAL600
! !INTERFACE:
subroutine AWRAL600_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use AWRAL600_lsmMod

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
!  The following is the list of variables written in the AWRAL600
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
    integer :: sr_ID
    integer :: sg_ID
    integer :: s0_ID
    integer :: ss_ID
    integer :: sd_ID
    integer :: mleaf_ID
    
    ! write the header of the restart file
    call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, &
                                       "AWRAL600", dim1=AWRAL600_struc(n)%nhru, &
				       dim2=AWRAL600_struc(n)%nhru, dimID=dimID, &
                                       output_format = trim(wformat))

    ! write the header for state variable sr
    call LIS_writeHeader_restart(ftn, n, dimID, sr_ID, "SR", &
                                 "volume of water in the surface water store", &
                                 "mm", vlevels=1, valid_min=0.0, valid_max=999.0)
    ! write the header for state variable sg
    call LIS_writeHeader_restart(ftn, n, dimID, sg_ID, "SG", &
                                 "groundwater storage in the unconfined aquifer", &
                                 "mm", vlevels=1, valid_min=0.0, valid_max=999.0)
    ! write the header for state variable s0
    call LIS_writeHeader_restart(ftn, n, dimID, s0_ID, "S0", &
                                 "water storage in the surface soil layer for each hru", &
                                 "mm", vlevels=AWRAL600_struc(n)%nhru , valid_min=0.0, valid_max=999.0, &
				 var_flag = "dim1")

 
    ! write the header for state variable ss
    call LIS_writeHeader_restart(ftn, n, dimID, ss_ID, "SS", &
                                 "water content of the shallow soil store for each hru", &
                                 "mm", vlevels=AWRAL600_struc(n)%nhru , valid_min=0.0, valid_max=999.0, &
				 var_flag = "dim1")
 
    ! write the header for state variable sd
    call LIS_writeHeader_restart(ftn, n, dimID, sd_ID, "SD", &
                                 "water content of the deep soil store for each hru", &
                                 "mm", vlevels=AWRAL600_struc(n)%nhru , valid_min=0.0, valid_max=999.0, &
				 var_flag = "dim1")
 
    ! write the header for state variable mleaf
    call LIS_writeHeader_restart(ftn, n, dimID, mleaf_ID, "MLEAF", &
                                 "leaf biomass", &
                                 "kg/m2", vlevels=AWRAL600_struc(n)%nhru , valid_min=0.0, valid_max=999.0, &
				 var_flag = "dim1")
 
    ! close header of restart file
    call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, AWRAL600_struc(n)%rstInterval)

    ! write state variables into restart file
    ! volume of water in the surface water store
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AWRAL600_struc(n)%awral600%sr, &
                              varid=sr_ID, dim=1, wformat=wformat)

    ! groundwater storage in the unconfined aquifer
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, AWRAL600_struc(n)%awral600%sg, &
                              varid=sg_ID, dim=1, wformat=wformat)

    ! water storage in the surface soil layer for each hru
    do l=1, AWRAL600_struc(n)%nhru
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = AWRAL600_struc(n)%awral600(t)%s0(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=s0_ID, dim=l, wformat=wformat)
    enddo
    ! water content of the shallow soil store for each hru
    do l=1, AWRAL600_struc(n)%nhru
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = AWRAL600_struc(n)%awral600(t)%ss(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=ss_ID, dim=l, wformat=wformat)
    enddo
    ! water content of the deep soil store for each hru
    do l=1, AWRAL600_struc(n)%nhru
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = AWRAL600_struc(n)%awral600(t)%sd(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=sd_ID, dim=l, wformat=wformat)
    enddo
    ! leaf biomass
    do l=1, AWRAL600_struc(n)%nhru
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = AWRAL600_struc(n)%awral600(t)%mleaf(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=mleaf_ID, dim=l, wformat=wformat)
    enddo
end subroutine AWRAL600_dump_restart
