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
! !ROUTINE: RUC37_writerst
! \label{RUC37_writerst}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   1/15/15: Shugong Wang; initial implementation for LIS 7 and RUC37
!
! !INTERFACE:
subroutine RUC37_writerst(n)
! !USES:
    use LIS_coreMod, only    : LIS_rc, LIS_masterproc
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber , LIS_verify
    use LIS_fileIOMod, only  : LIS_create_output_directory, &
                               LIS_create_restart_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN
    use RUC37_lsmMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: n
!
! !DESCRIPTION:
!  This program writes restart files for RUC37.
!  This includes all relevant water/energy storage and tile information.
!
!  The routines invoked are:
! \begin{description}
! \item[LIS\_create\_output\_directory](\ref{LIS_create_output_directory}) \newline
!  creates a timestamped directory for the restart files
! \item[LIS\_create\_restart\_filename](\ref{LIS_create_restart_filename}) \newline
!  generates a timestamped restart filename
! \item[RUC37\_dump\_restart](\ref{RUC37_dump_restart}) \newline
!   writes the RUC37 variables into the restart file
! \end{description}
!EOP

    character(len=LIS_CONST_PATH_LEN) :: filen
    character*20  :: wformat
    logical       :: alarmCheck
    integer       :: ftn
    integer       :: status
    
    ! set restart alarm
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "RUC37 restart alarm")
    
    ! set restart file format (read from LIS configration file_
    wformat = trim(RUC37_struc(n)%rformat)
    
    if(alarmCheck .or. (LIS_rc%endtime ==1)) then
        If (LIS_masterproc) Then
            call LIS_create_output_directory("SURFACEMODEL")
            call LIS_create_restart_filename(n, filen, "SURFACEMODEL", "RUC37",&
                                             wformat=wformat)
            if(wformat .eq. "binary") then
                ftn = LIS_getNextUnitNumber()
                open(ftn,file=filen,status="unknown", form="unformatted")
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF4)
                status = nf90_create(path=filen, cmode=nf90_hdf5, ncid = ftn)
                call LIS_verify(status,"Error in nf90_open in RUC37_writerst")
#endif
#if (defined USE_NETCDF3)
                status = nf90_create(Path = filen, cmode = nf90_clobber, ncid = ftn)
                call LIS_verify(status, "Error in nf90_open in RUC37_writerst")
#endif
             endif
        endif
    
        call RUC37_dump_restart(n, ftn, wformat)
    
        if (LIS_masterproc) then
            if(wformat .eq. "binary") then
                call LIS_releaseUnitNumber(ftn)
            elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
                status = nf90_close(ftn)
                call LIS_verify(status, "Error in nf90_close in RUC37_writerst")
#endif
            endif
            write(LIS_logunit, *) "RUC37 archive restart written: ", trim(filen)
        endif
    endif
end subroutine RUC37_writerst

!BOP
!
! !ROUTINE: RUC37_dump_restart
! \label{RUC37_dump_restart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  1/15/15: Shugong Wang, initial implementation for LIS 7 and RUC37
! !INTERFACE:
subroutine RUC37_dump_restart(n, ftn, wformat)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use RUC37_lsmMod

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
!  The following is the list of variables written in the RUC37
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
!    smfr                       - RUC37 soil ice content (m3 m-3) [m^3 m-3]
!    keepfr                     - RUC37 frozen soil flag: 0 or 1
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
    integer :: emiss_ID
    integer :: ch_ID
    integer :: cm_ID
    integer :: sneqv_ID
    integer :: snowh_ID
    integer :: snowc_ID
    integer :: canwat_ID
    integer :: alb_ID
    integer :: smc_ID
    integer :: sho_ID
    integer :: stc_ID
    integer :: smfr_ID
    integer :: keepfr_ID
    integer :: tskin_ID
    integer :: qvg_ID
    integer :: qsfc_ID
    integer :: qcg_ID
    integer :: qsg_ID
    integer :: snt75cm_ID
    integer :: tsnav_ID
    integer :: soilm_ID
    integer :: smroot_ID
    
    ! write the header of the restart file
    call LIS_writeGlobalHeader_restart(ftn, n, LIS_rc%lsm_index, &
                                       "RUC37", dim1=RUC37_struc(n)%nsoil, &
                                       dim2=4, dimID=dimID, &
                                       output_format = trim(wformat))

    ! write the header for state variable emiss
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, emiss_ID, "EMISS", &
                                 "surface emissivity (0.0 - 1.0).", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable ch
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, ch_ID, "CH", &
                                 "exchange coefficient for head and moisture (m s-1).", &
                                 "s m-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable cm
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, cm_ID, "CM", &
                                 "exchange coefficient for momentum (m s-1).", &
                                 "s m-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable sneqv
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sneqv_ID, "SNEQV", &
                                 "water equivalent of accumulated snow depth (m).", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable snowh
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, snowh_ID, "SNOWH", &
                                 "physical snow depth (m).", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable snowc
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, snowc_ID, "SNOWC", &
                                 "fractional snow cover ( fraction [0.0-1.0] )", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable canwat
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, canwat_ID, "CANWAT", &
                                 "canopy moisture content (kg m-2)", &
                                 "kg m-2", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable alb
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, alb_ID, "ALB", &
                                 "surface albedo including possible snow-cover effect.  this is set in lsmruc,", &
                                 "-", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable smc
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, smc_ID, "SMC", &
                                 "total soil moisture content (m3 m-3)", &
                                 "m^3 m-3", vlevels=RUC37_struc(n)%nsoil, valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")
 
    ! write the header for state variable sho
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, sho_ID, "SHO", &
                                 "liquid soil moisture content (m3 m-3)", &
                                 "m^3 m-3", vlevels=RUC37_struc(n)%nsoil , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1") 
 
    ! write the header for state variable stc
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, stc_ID, "STC", &
                                 "soil temperature (k)", &
                                 "K", vlevels=RUC37_struc(n)%nsoil , valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")

    ! write the header for state variable smfr
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and
    !valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, smfr_ID, "SMFR", &
                                 "soil ice content (m3 m-3)", &
                                 "m^3 m-3", vlevels=RUC37_struc(n)%nsoil, valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")

    ! write the header for state variable keepfr
    !TODO: check dimension of the state variable following "vlevels="
    !TODO: replace -99999 and 99999 with correct values for valid_min and
    !valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, keepfr_ID, "KEEPFR", &
                                 "frozen soil flag )", &
                                 " ", vlevels=RUC37_struc(n)%nsoil, valid_min=-99999.0, valid_max=99999.0, &
                                 var_flag = "dim1")
 
    ! write the header for state variable tskin
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, tskin_ID, "TSKIN", &
                                 "skin temperature (k)", &
                                 "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable qvg
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qvg_ID, "QVG", &
                                 "mixing ratio at the surface ( kg kg{-1} )", &
                                 "kg kg-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable qsfc
    !TODO: replace -99999 and 99999 with correct values for valid_min and
    !valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qsfc_ID, "QSFC", &
                                 "specific humidity at the surface ( kg kg{-1} )", &
                                 "kg kg-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)

    ! write the header for state variable qcg
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qcg_ID, "QCG", &
                                 "cloud water mixing ratio at the surface ( kg kg{-1} )", &
                                 "kg kg-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable qsg
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, qsg_ID, "QSG", &
                                 "surface water vapor mixing ratio at satration (kg kg-1)", &
                                 "kg kg-1", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable snt75cm
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, snt75cm_ID, "SNT75CM", &
                                 "snow temperature at 7.5 cm depth (k)", &
                                 "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable tsnav
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, tsnav_ID, "TSNAV", &
                                 "average snow temperature in k", &
                                 "K", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable soilm
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, soilm_ID, "SOILM", &
                                 "total soil column moisture content, frozen and unfrozen ( m )", &
                                 "m", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! write the header for state variable smroot
    !TODO: replace -99999 and 99999 with correct values for valid_min and valid_max
    call LIS_writeHeader_restart(ftn, n, dimID, smroot_ID, "SMROOT", &
                                 "available soil moisture in the root zone ( fraction [smcwlt-smcmax]", &
                                 "m^3 m-3", vlevels=1, valid_min=-99999.0, valid_max=99999.0)
    ! close header of restart file
    call LIS_closeHeader_restart(ftn, n, LIS_rc%lsm_index, dimID, RUC37_struc(n)%rstInterval)

    ! write state variables into restart file
    ! surface emissivity (0.0 - 1.0).
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%emiss, &
                              varid=emiss_ID, dim=1, wformat=wformat)

    ! exchange coefficient for head and moisture (m s-1).
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%ch, &
                              varid=ch_ID, dim=1, wformat=wformat)

    ! exchange coefficient for momentum (m s-1).
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%cm, &
                              varid=cm_ID, dim=1, wformat=wformat)

    ! water equivalent of accumulated snow depth (m).
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%sneqv, &
                              varid=sneqv_ID, dim=1, wformat=wformat)

    ! physical snow depth (m).
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%snowh, &
                              varid=snowh_ID, dim=1, wformat=wformat)

    ! fractional snow cover ( fraction [0.0-1.0] )
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%snowc, &
                              varid=snowc_ID, dim=1, wformat=wformat)

    ! canopy moisture content (kg m-2)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%canwat, &
                              varid=canwat_ID, dim=1, wformat=wformat)

    ! surface albedo including possible snow-cover effect.  this is set in lsmruc,
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%alb, &
                              varid=alb_ID, dim=1, wformat=wformat)

    ! total soil moisture content (m3 m-3)
    do l=1, RUC37_struc(n)%nsoil   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = RUC37_struc(n)%ruc37(t)%smc(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=smc_ID, dim=l, wformat=wformat)
    enddo
    ! liquid soil moisture content (m3 m-3)
    do l=1, RUC37_struc(n)%nsoil   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = RUC37_struc(n)%ruc37(t)%sho(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=sho_ID, dim=l, wformat=wformat)
    enddo
    ! soil temperature (k)
    do l=1, RUC37_struc(n)%nsoil   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = RUC37_struc(n)%ruc37(t)%stc(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=stc_ID, dim=l, wformat=wformat)
    enddo
    ! soil ice content (m3 m-3)
    do l=1, RUC37_struc(n)%nsoil   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = RUC37_struc(n)%ruc37(t)%smfr(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=smfr_ID, dim=l, wformat=wformat)
    enddo
    ! frozen soil flag (m3 m-3)
    do l=1, RUC37_struc(n)%nsoil   ! TODO: check loop
        tmptilen = 0
        do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            tmptilen(t) = RUC37_struc(n)%ruc37(t)%keepfr(l)
        enddo
        call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, &
                                  varid=keepfr_ID, dim=l, wformat=wformat)
    enddo

    ! skin temperature (k)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%tskin, &
                              varid=tskin_ID, dim=1, wformat=wformat)

    ! mixing ratio at the surface ( kg kg{-1} )
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%qvg, &
                              varid=qvg_ID, dim=1, wformat=wformat)
    ! specific humidity at the surface ( kg kg{-1} )
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%qsfc, &
                              varid=qsfc_ID, dim=1, wformat=wformat)
    ! effective cloud water mixing ratio at the surface ( kg kg{-1} )
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%qcg, &
                              varid=qcg_ID, dim=1, wformat=wformat)

    ! surface water vapor mixing ratio at satration (kg kg-1)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%qsg, &
                              varid=qsg_ID, dim=1, wformat=wformat)

    ! snow temperature at 7.5 cm depth (k)
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%snt75cm, &
                              varid=snt75cm_ID, dim=1, wformat=wformat)

    ! average snow temperature in k
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%tsnav, &
                              varid=tsnav_ID, dim=1, wformat=wformat)

    ! total soil column moisture content, frozen and unfrozen ( m )
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%soilm, &
                              varid=soilm_ID, dim=1, wformat=wformat)

    ! available soil moisture in the root zone ( fraction [smcwlt-smcmax]
    call LIS_writevar_restart(ftn, n, LIS_rc%lsm_index, RUC37_struc(n)%ruc37%smroot, &
                              varid=smroot_ID, dim=1, wformat=wformat)

end subroutine RUC37_dump_restart
