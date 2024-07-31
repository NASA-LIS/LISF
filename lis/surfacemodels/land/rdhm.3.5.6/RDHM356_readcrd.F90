!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: RDHM356_readcrd
! \label{RDHM356_readcrd}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   11/5/13 : Shugong Wang, initial implementation for LIS 7 and RDHM356
!
! !INTERFACE:
subroutine RDHM356_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit , LIS_verify
    use RDHM356_lsmMod, only : RDHM356_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to RDHM356 model from
!  the LIS configuration file.
!
!EOP
    implicit none

    integer      :: rc 
    integer      :: n, i
    character*10 :: time 
    character*6  :: str_i
    real         :: init_ratio

    write(LIS_logunit, *) "Start reading LIS configuration file for RDHM356 model"
    
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "RDHM356 model timestep: not defined")
        call LIS_parseTimeString(time, RDHM356_struc(n)%ts)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc,"RDHM356 restart output interval: not defined")
        call LIS_parseTimeString(time, RDHM356_struc(n)%rstInterval)
    enddo
    
    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number of desired soil layers for total and liquid soil moisture
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 NDINTW:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%NDINTW, rc=rc)
        call LIS_verify(rc, "RDHM356 NDINTW: not defined")
    enddo
 
    do n=1, LIS_rc%nnest
        allocate(RDHM356_struc(n)%DSINTW(RDHM356_struc(n)%NDINTW))
    enddo

    ! number of desired soil layers for soil temperature
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 NDSINT:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%NDSINT, rc=rc)
        call LIS_verify(rc, "RDHM356 NDSINT: not defined")
    enddo
 
    do n=1, LIS_rc%nnest
        allocate(RDHM356_struc(n)%DSINT(RDHM356_struc(n)%NDSINT))
    enddo

    ! observation height of temperature of humidity
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 TempHeight:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%TempHeight, rc=rc)
        call LIS_verify(rc, "RDHM356 TempHeight: not defined")
    enddo
 
    ! observation height of wind
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 WindHeight:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%WindHeight, rc=rc)
        call LIS_verify(rc, "RDHM356 WindHeight: not defined")
    enddo
 
    ! simulation time interval of SAC model and Snow-17
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 DT_SAC_SNOW17:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%DT_SAC_SNOW17, rc=rc)
        call LIS_verify(rc, "RDHM356 DT_SAC_SNOW17: not defined")
    enddo
 
    ! simulation time interval of frozen soil model
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 DT_FRZ:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%DT_FRZ, rc=rc)
        call LIS_verify(rc, "RDHM356 DT_FRZ: not defined")
    enddo
 
    ! version number of frozen soil model. 1: old version, 2: new version
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 FRZ_VER_OPT:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%FRZ_VER_OPT, rc=rc)
        call LIS_verify(rc, "RDHM356 FRZ_VER_OPT: not defined")
    enddo

    ! option for SACHTET. If SACHTET_OPT=1, run SACHTET, otherwise, don't run
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 SACHTET_OPT:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%SACHTET_OPT, rc=rc)
        call LIS_verify(rc, "RDHM356 SACHTET_OPT: not defined")
    enddo
 
 
    ! option for snow-17. If SNOW17_OPT=1, use SNOW-17, otherwise, don't use
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 SNOW17_OPT:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%SNOW17_OPT, rc=rc)
        call LIS_verify(rc, "RDHM356 SNOW17_OPT: not defined")
    enddo
 
    ! number of soil types
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 NSTYP:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%NSTYP, rc=rc)
        call LIS_verify(rc, "RDHM356 NSTYP: not defined")
    enddo
 
    ! number of vegetation types
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 NVTYP:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%NVTYP, rc=rc)
        call LIS_verify(rc, "RDHM356 NVTYP: not defined")
    enddo
 
    ! normalization flag for total and liquid soil moisture output (1-normalized, 0-not)
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 NORMALIZE:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%NORMALIZE, rc=rc)
        call LIS_verify(rc, "RDHM356 NORMALIZE: not defined")
    enddo
 
    ! thickness of desired soil layers for liquid and total soil moisture
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 DSINTW:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, RDHM356_struc(n)%NDINTW
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%DSINTW(i), rc=rc)
            call LIS_verify(rc, 'RDHM356 DSINTW: not defined')
        enddo
    enddo
 
    ! thickness of desired soil layers for soil temperature
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 DSINT:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, RDHM356_struc(n)%NDSINT
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%DSINT(i), rc=rc)
            call LIS_verify(rc, 'RDHM356 DSINT: not defined')
        enddo
    enddo
 
    ! adjustment of PET for 12 months
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 PETADJ_MON:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, 12
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%PETADJ_MON(i), rc=rc)
            call LIS_verify(rc, 'RDHM356 PETADJ_MON: not defined')
        enddo
    enddo
 
    ! default=0.12 Zilitinkevich
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 CZIL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%CZIL, rc=rc)
        call LIS_verify(rc, "RDHM356 CZIL: not defined")
    enddo
 
    ! FXEXP(fxexp),(default=2.0) bare soil
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 FXEXP:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%FXEXP, rc=rc)
        call LIS_verify(rc, "RDHM356 FXEXP: not defined")
    enddo
 
    ! RCMAX,(default=5000s/m) maximum stomatal resistance
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 vegRCMAX:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%vegRCMAX, rc=rc)
        call LIS_verify(rc, "RDHM356 vegRCMAX: not defined")
    enddo
 
    ! TOPT,(default=298K)optimum air
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 TOPT:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%TOPT, rc=rc)
        call LIS_verify(rc, "RDHM356 TOPT: not defined")
    enddo
 
    ! plant coef. default pc = -1, 0.6 - 0.8
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 PC:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%PC, rc=rc)
        call LIS_verify(rc, "RDHM356 PC: not defined")
    enddo
 
    ! if PET_OPT = 0, use non Penmann-based ETP;if penpt > 0 empirical Penmann equation; if penpt < 0, use energy based Pennman
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 PET_OPT:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%PET_OPT, rc=rc)
        call LIS_verify(rc, "RDHM356 PET_OPT: not defined")
    enddo
 
    ! default=1 means noah option,this constant allows selection of tension water redistribution option, 
    ! if rdst = 0 (ohd), use OHD version of SRT subroutine this SRT uses reference gradient instead an actual.
    ! if rdst = 1 ( noah), use Noah version of SRT subroutine
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 RDST:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%RDST, rc=rc)
        call LIS_verify(rc, "RDHM356 RDST: not defined")
    enddo
 
    ! this constant allows change of RCMIN (0.5)
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 thresholdRCMIN:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%thresholdRCMIN, rc=rc)
        call LIS_verify(rc, "RDHM356 thresholdRCMIN: not defined")
    enddo
 
    ! reference wind speed for PET adjustment (4 m s-1)
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 SFCREF:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%SFCREF, rc=rc)
        call LIS_verify(rc, "RDHM356 SFCREF: not defined")
    enddo
 
    ! Ek-Chen evaporation threshold switch. Bare soil evaporation option changes according to greenness.
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 BAREADJ:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%BAREADJ, rc=rc)
        call LIS_verify(rc, "RDHM356 BAREADJ: not defined")
    enddo
 
    ! switch variable change liquid water freezing version, 0: Victor's version, 1: Eric's version
    call ESMF_ConfigFindLabel(LIS_config, "RDHM356 SNOW17_SWITCH:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%SNOW17_SWITCH, rc=rc)
        call LIS_verify(rc, "RDHM356 SNOW17_SWITCH: not defined")
    enddo

    !The following lines hard code the LDT NetCDF variable names. 
    do n=1, LIS_rc%nnest
        RDHM356_struc(n)%LDT_ncvar_Tair_min  = 'RDHM356_TAIR_MIN'
        RDHM356_struc(n)%LDT_ncvar_Tair_max  = 'RDHM356_TAIR_MAX'
        !RDHM356_struc(n)%LDT_ncvar_PET_MON   = 'RDHM356_PET_MON'
        RDHM356_struc(n)%LDT_ncvar_PET_MON   = 'PET'
        !RDHM356_struc(n)%LDT_ncvar_GRN_MON   = 'RDHM356_GRN_MON'
        RDHM356_struc(n)%LDT_ncvar_GRN_MON   = 'GREENNESS'
        RDHM356_struc(n)%LDT_ncvar_SoilAlb   = 'RDHM356_SOILALB'
        RDHM356_struc(n)%LDT_ncvar_SnowAlb   = 'MXSNALBEDO'
        RDHM356_struc(n)%LDT_ncvar_SOILTYP   = 'RDHM356_SOILTEXT'
        RDHM356_struc(n)%LDT_ncvar_VEGETYP   = 'RDHM356_VEGETYP'
        RDHM356_struc(n)%LDT_ncvar_UZTWM     = 'RDHM356_UZTWM'
        RDHM356_struc(n)%LDT_ncvar_UZFWM     = 'RDHM356_UZFWM'
        RDHM356_struc(n)%LDT_ncvar_UZK       = 'RDHM356_UZK'
        RDHM356_struc(n)%LDT_ncvar_PCTIM     = 'RDHM356_PCTIM'
        RDHM356_struc(n)%LDT_ncvar_ADIMP     = 'RDHM356_ADIMP'
        RDHM356_struc(n)%LDT_ncvar_RIVA      = 'RDHM356_RIVA'
        RDHM356_struc(n)%LDT_ncvar_ZPERC     = 'RDHM356_ZPERC'
        RDHM356_struc(n)%LDT_ncvar_REXP      = 'RDHM356_REXP'
        RDHM356_struc(n)%LDT_ncvar_LZTWM     = 'RDHM356_LZTWM'
        RDHM356_struc(n)%LDT_ncvar_LZFSM     = 'RDHM356_LZFSM'
        RDHM356_struc(n)%LDT_ncvar_LZFPM     = 'RDHM356_LZFPM'
        RDHM356_struc(n)%LDT_ncvar_LZSK      = 'RDHM356_LZSK'
        RDHM356_struc(n)%LDT_ncvar_LZPK      = 'RDHM356_LZPK'
        RDHM356_struc(n)%LDT_ncvar_PFREE     = 'RDHM356_PFREE'
        RDHM356_struc(n)%LDT_ncvar_SIDE      = 'RDHM356_SIDE'
        RDHM356_struc(n)%LDT_ncvar_RSERV     = 'RDHM356_RSERV'
        RDHM356_struc(n)%LDT_ncvar_EFC       = 'RDHM356_EFC'
        RDHM356_struc(n)%LDT_ncvar_TBOT      = 'TBOT'
        RDHM356_struc(n)%LDT_ncvar_RSMAX     = 'RDHM356_RSMAX'
        RDHM356_struc(n)%LDT_ncvar_CKSL      = 'RDHM356_CKSL'
        RDHM356_struc(n)%LDT_ncvar_ZBOT      = 'RDHM356_ZBOT'
        RDHM356_struc(n)%LDT_ncvar_vegRCMIN  = 'RDHM356_RCMIN'
        RDHM356_struc(n)%LDT_ncvar_climRCMIN = 'RDHM356_RCMINCLIM'
        RDHM356_struc(n)%LDT_ncvar_RGL       = 'RDHM356_RGL'
        RDHM356_struc(n)%LDT_ncvar_HS        = 'RDHM356_HS'
        RDHM356_struc(n)%LDT_ncvar_LAI       = 'RDHM356_MAXLAI'
        RDHM356_struc(n)%LDT_ncvar_D50       = 'RDHM356_D50'
        RDHM356_struc(n)%LDT_ncvar_CROOT     = 'RDHM356_CROOT'
        RDHM356_struc(n)%LDT_ncvar_Z0        = 'RDHM356_Z0'
        RDHM356_struc(n)%LDT_ncvar_CLAY      = 'RDHM356_CLAY'
        RDHM356_struc(n)%LDT_ncvar_SAND      = 'RDHM356_SAND'
        RDHM356_struc(n)%LDT_ncvar_SATDK     = 'RDHM356_SATDK'
        RDHM356_struc(n)%LDT_ncvar_ALON      = 'lon'
        RDHM356_struc(n)%LDT_ncvar_ALAT      = 'RDHM356_ALAT'
        RDHM356_struc(n)%LDT_ncvar_SCF       = 'RDHM356_SCF'
        RDHM356_struc(n)%LDT_ncvar_MFMAX     = 'RDHM356_MFMAX'
        RDHM356_struc(n)%LDT_ncvar_MFMIN     = 'RDHM356_MFMIN'
        RDHM356_struc(n)%LDT_ncvar_NMF       = 'RDHM356_NMF'
        RDHM356_struc(n)%LDT_ncvar_UADJ      = 'RDHM356_UADJ'
        RDHM356_struc(n)%LDT_ncvar_SI        = 'RDHM356_SI'
        RDHM356_struc(n)%LDT_ncvar_MBASE     = 'RDHM356_MBASE'
        RDHM356_struc(n)%LDT_ncvar_PXTEMP    = 'RDHM356_PXTEMP'
        RDHM356_struc(n)%LDT_ncvar_PLWHC     = 'RDHM356_PLWHC'
        RDHM356_struc(n)%LDT_ncvar_TIPM      = 'RDHM356_TIPM'
        RDHM356_struc(n)%LDT_ncvar_GM        = 'RDHM356_PGM'
        RDHM356_struc(n)%LDT_ncvar_ELEV      = 'RDHM356_ELEV'
        RDHM356_struc(n)%LDT_ncvar_LAEC      = 'RDHM356_LAEC'
        RDHM356_struc(n)%LDT_ncvar_ADC       = 'ADC'
    enddo
    ! restart run, read restart file
    do n=1,LIS_rc%nnest
        RDHM356_struc(n)%rformat='netcdf' ! default restart format  
    enddo
        Call ESMF_ConfigFindLabel(LIS_config, "RDHM356 tmxmn directory:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%tmxmn_dir, rc=rc)
            call LIS_verify(rc, "RDHM356 tmxmn directory: not defined")
        enddo
    if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, "RDHM356 restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, "RDHM356 restart file: not defined")
        enddo
        
        Call ESMF_ConfigFindLabel(LIS_config, "RDHM356 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, "RDHM356 restart file format: not defined")
        enddo
    ! cold start run, read initial state variables
    else 
        ! upper zone tension water storage content
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial UZTWC (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_UZTWC_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial UZTWC: not defined")
        enddo

        ! upper zone free water storage content
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial UZFWC (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_UZFWC_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial UZFWC: not defined")
        enddo

        ! lower zone tension water storage content
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial LZTWC (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_LZTWC_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial LZTWC: not defined")
        enddo

        ! lower zone primary free water storage content
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial LZFPC (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_LZFPC_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial LZFPC: not defined")
        enddo

        ! lower zone supplemental free water storage content
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial LZFSC (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_LZFSC_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial LZFSC: not defined")
        enddo

        ! additional impervious area content
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial ADIMC (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_ADIMC_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial ADIMC: not defined")
        enddo

        ! first soil layer temperature
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial TS0:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_TS0, rc=rc)
            call LIS_verify(rc, "RDHM356 initial TS0: not defined")
        enddo

!        ! second soil layer temperature
!        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial TS1:", rc = rc)
!        do n=1,LIS_rc%nnest
!            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_TS1, rc=rc)
!            call LIS_verify(rc, "RDHM356 initial TS1: not defined")
!        enddo
!
!        ! third soil layer temperature
!        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial TS2:", rc = rc)
!        do n=1,LIS_rc%nnest
!            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_TS2, rc=rc)
!            call LIS_verify(rc, "RDHM356 initial TS2: not defined")
!        enddo
!
!        ! fourth soil layer temperature
!        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial TS3:", rc = rc)
!        do n=1,LIS_rc%nnest
!            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_TS3, rc=rc)
!            call LIS_verify(rc, "RDHM356 initial TS3: not defined")
!        enddo
!
!        ! fifth soil layer temperature
!        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial TS4:", rc = rc)
!        do n=1,LIS_rc%nnest
!            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_TS4, rc=rc)
!            call LIS_verify(rc, "RDHM356 initial TS4: not defined")
!        enddo

        ! unfrozen upper zone tension water
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial UZTWH (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_UZTWH_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial UZTWH: not defined")
        enddo

        ! unfrozen uppeer zone free water
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial UZFWH (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_UZFWH_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial UZFWH: not defined")
        enddo

        ! unfrozen lower zone tension water
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial LZTWH (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_LZTWH_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial LZTWH: not defined")
        enddo

        ! unfrozen lower zone supplemental free water
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial LZFSH (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_LZFSH_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial LZFSH: not defined")
        enddo

        ! unfrozen lower zone primary free water
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial LZFPH (ratio):", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_LZFPH_ratio, rc=rc)
            call LIS_verify(rc, "RDHM356 initial LZFPH: not defined")
        enddo

        ! volumetric content of total soil moisture at each layer
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial SMC:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, 6 !!! TODO: check dimension here
                call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_SMC(i), rc=rc)
            end do
            call LIS_verify(rc, "RDHM356 initial SMC: not defined")
        enddo

        ! volumetric content of liquid soil moisture at each layer
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial SH2O:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, 6 !!! TODO: check dimension here
                call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_SH2O(i), rc=rc)
            end do
            call LIS_verify(rc, "RDHM356 initial SH2O: not defined")
        enddo

        ! snow water equivalent without liquid water
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial WE:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_WE, rc=rc)
            call LIS_verify(rc, "RDHM356 initial WE: not defined")
        enddo

        ! liquid water in snow
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial LIQW:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_LIQW, rc=rc)
            call LIS_verify(rc, "RDHM356 initial LIQW: not defined")
        enddo

        ! negative snow heat
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial NEGHS:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_NEGHS, rc=rc)
            call LIS_verify(rc, "RDHM356 initial NEGHS: not defined")
        enddo

        ! antecedent temperature index
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial TINDEX:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_TINDEX, rc=rc)
            call LIS_verify(rc, "RDHM356 initial TINDEX: not defined")
        enddo

        ! cumulated snow water including liquid
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial ACCMAX:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_ACCMAX, rc=rc)
            call LIS_verify(rc, "RDHM356 initial ACCMAX: not defined")
        enddo

        ! snow depth
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial SNDPT:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_SNDPT, rc=rc)
            call LIS_verify(rc, "RDHM356 initial SNDPT: not defined")
        enddo

        ! average snow temperature
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial SNTMP:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_SNTMP, rc=rc)
            call LIS_verify(rc, "RDHM356 initial SNTMP: not defined")
        enddo

        ! the last highest snow water equivalent before any snow fall
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial SB:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_SB, rc=rc)
            call LIS_verify(rc, "RDHM356 initial SB: not defined")
        enddo

        ! internal snow state during melt & new snow fall (checked with Victor)
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial SBAESC:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_SBAESC, rc=rc)
            call LIS_verify(rc, "RDHM356 initial SBAESC: not defined")
        enddo

        ! internal snow state during melt & new snow fall (checked with Victor)
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial SBWS:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_SBWS, rc=rc)
            call LIS_verify(rc, "RDHM356 initial SBWS: not defined")
        enddo

        ! snow liquid water attenuation storage
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial STORAGE:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_STORAGE, rc=rc)
            call LIS_verify(rc, "RDHM356 initial STORAGE: not defined")
        enddo

        ! adjusted areal snow cover fraction
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial AEADJ:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_AEADJ, rc=rc)
            call LIS_verify(rc, "RDHM356 initial AEADJ: not defined")
        enddo

        ! array of lagged liquid water values
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial EXLAG:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, 7 !!! TODO: check dimension here
                call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_EXLAG(i), rc=rc)
            end do
            call LIS_verify(rc, "RDHM356 initial EXLAG: not defined")
        enddo

        ! number of ordinates in lagged liquid water array (EXLAG)
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial NEXLAG:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_NEXLAG, rc=rc)
            call LIS_verify(rc, "RDHM356 initial NEXLAG: not defined")
        enddo

        ! air temperature of previous time step
        call ESMF_ConfigFindLabel(LIS_config, "RDHM356 initial TA_PREV:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RDHM356_struc(n)%init_TA_PREV, rc=rc)
            call LIS_verify(rc, "RDHM356 initial TA_PREV: not defined")
        enddo

    end if
     
    write(LIS_logunit, *) "Finish reading LIS configuration file for RDHM356 model"
     
end subroutine RDHM356_readcrd

subroutine check_init_ratio(init_ratio, label)
    use LIS_logMod, only     : LIS_logunit 
    implicit none
    real, intent(in) :: init_ratio
    character(len=*) :: label

    if(init_ratio<0.0 .or. init_ratio>1.0) then
        write(LIS_logunit, *) "RDHM356 initialization failed, check the ragne of "
        write(LIS_logunit, *) label
        write(LIS_logunit, *) "input value :", init_ratio, " out of (0.0 - 1.0)"
    endif
end subroutine check_init_ratio
