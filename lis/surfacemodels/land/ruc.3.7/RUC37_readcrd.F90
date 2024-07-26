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
! !ROUTINE: RUC37_readcrd
! \label{RUC37\_readcrd}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   1/15/15 : Shugong Wang, initial implementation for LIS 7 and RUC37
!
! !INTERFACE:
subroutine RUC37_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit , LIS_verify
    use RUC37_lsmMod, only   : RUC37_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to RUC37 model from
!  the LIS configuration file.
!
!EOP
    implicit none

    integer      :: rc 
    integer      :: n, i
    character*10 :: time 
    character*6  :: str_i
 
    write(LIS_logunit, *) "Start reading LIS configuration file for RUC37 model"
    
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "RUC37 model timestep: not defined")
        call LIS_parseTimeString(time, RUC37_struc(n)%ts)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc,"RUC37 restart output interval: not defined")
        call LIS_parseTimeString(time, RUC37_struc(n)%rstInterval)
    enddo
    
    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number of soil levels.
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 number of soil levels:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%nsoil, rc=rc)
        call LIS_verify(rc, "RUC37 nsoil: not defined")
    enddo
 
    ! allocate memory for soil_layer_thickness using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(RUC37_struc(n)%soil_layer_thickness(RUC37_struc(n)%nsoil))
    enddo
    ! allocate memory for init_smc using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(RUC37_struc(n)%init_smc(RUC37_struc(n)%nsoil))
    enddo
    ! allocate memory for init_sho using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(RUC37_struc(n)%init_sho(RUC37_struc(n)%nsoil))
    enddo
    ! allocate memory for init_stc using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(RUC37_struc(n)%init_stc(RUC37_struc(n)%nsoil))
    enddo
    ! allocate memory for init_smfr using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(RUC37_struc(n)%init_smfr(RUC37_struc(n)%nsoil))
    enddo
    ! allocate memory for init_keepfr using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(RUC37_struc(n)%init_keepfr(RUC37_struc(n)%nsoil))
    enddo
 
    ! time step (seconds).
    !call ESMF_ConfigFindLabel(LIS_config, "RUC37 dt:", rc = rc)
    do n=1, LIS_rc%nnest
        !call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%dt, rc=rc)
        !call LIS_verify(rc, "RUC37 dt: not defined")
         RUC37_struc(n)%dt =  RUC37_struc(n)%ts 
    enddo
 
    ! thicknesses of each soil level (m)
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 soil level depth:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, RUC37_struc(n)%nsoil
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%soil_layer_thickness(i), rc=rc)
            call LIS_verify(rc, 'RUC37 soil_layer_thickness: not defined')
        enddo
    enddo
 
    ! .true. to use table values for albbck, shdfac, and z0brd; .false. to use values for albbck, shdfac, and z0brd as set in this driver routine
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 use local parameters:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%use_local_param, rc=rc)
        call LIS_verify(rc, "RUC37 use local parameters: not defined")
    enddo
 
    ! if rdlai2d == .true., then the xlai value that we pass to lsmruc will be used. if rdlai2d == .false., then xlai will be computed within lsmruc, from table minimum and maximum values in vegparm.tbl, and the current green vegetation fraction.
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 use 2D LAI map:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%use_2d_lai_map, rc=rc)
        call LIS_verify(rc, "RUC37 use 2D LAI map: not defined")
    enddo
 
    ! if usemonalb == .true., then the alb value passed to lsmruc will be used as the background snow-free albedo term.  if usemonalb == .false., then alb will be computed within lsmruc from minimum and maximum values in vegparm.tbl, and the current green vegetation fraction.
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 use monthly albedo map:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%use_monthly_albedo_map, rc=rc)
        call LIS_verify(rc, "RUC37 use monthly albedo map: not defined")
    enddo
 
    ! option to turn on (iz0tlnd=1) or off (iz0tlnd=0) the vegetation-category-dependent calculation of the zilitinkivich coefficient czil in the sfcdif subroutines.
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 option_iz0tlnd:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%option_iz0tlnd, rc=rc)
        call LIS_verify(rc, "RUC37 option_iz0tlnd: not defined")
    enddo
 
    ! option to use previous (sfcdif_option=0) or updated (sfcdif_option=1) version of sfcdif subroutine.
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 option_sfcdif:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%option_sfcdif, rc=rc)
        call LIS_verify(rc, "RUC37 option_sfcdif: not defined")
    enddo
 
    ! noah model landuse parameter table
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 landuse_tbl_name:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%landuse_tbl_name, rc=rc)
        call LIS_verify(rc, "RUC37 landuse_tbl_name: not defined")
    enddo
 
    ! noah model soil parameter table
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 soil_tbl_name:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%soil_tbl_name, rc=rc)
        call LIS_verify(rc, "RUC37 soil_tbl_name: not defined")
    enddo
    
    ! noah model soil parameter table
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 gen_tbl_name:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%gen_tbl_name, rc=rc)
        call LIS_verify(rc, "RUC37 gen_tbl_name: not defined")
    enddo
 
    ! landuse classification scheme
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 landuse_scheme_name:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%landuse_scheme_name, rc=rc)
        call LIS_verify(rc, "RUC37 landuse_scheme_name: not defined")
    enddo
 
    ! soil classification scheme
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 soil_scheme_name:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%soil_scheme_name, rc=rc)
        call LIS_verify(rc, "RUC37 soil_scheme_name: not defined")
    enddo
 
    ! number of water category in llanduse classification
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 water_class_num:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%water_class_num, rc=rc)
        call LIS_verify(rc, "RUC37 water_class_num: not defined")
    enddo
 
    ! number of ice category in llanduse classification
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 ice_class_num:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%ice_class_num, rc=rc)
        call LIS_verify(rc, "RUC37 ice_class_num: not defined")
    enddo
 
    ! number of urban category in llanduse classification
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 urban_class_num:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%urban_class_num, rc=rc)
        call LIS_verify(rc, "RUC37 urban_class_num: not defined")
    enddo

    ! reference height of temperature and humidity
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 zlvl:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%zlvl, rc=rc)
        call LIS_verify(rc, "RUC37 zlvl: not defined")
    enddo

    ! reference height of wind
    call ESMF_ConfigFindLabel(LIS_config, "RUC37 zlvl_wind:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%zlvl_wind, rc=rc)
        call LIS_verify(rc, "RUC37 zlvl_wind: not defined")
    enddo
 

    ! The following lines hard code the LDT NetCDF variable names. 
    do n=1, LIS_rc%nnest
!        RUC37_struc(n)%LDT_ncvar_vegetype    = 'RUC37_VEGETYPE'
!        RUC37_struc(n)%LDT_ncvar_soiltype    = 'RUC37_SOILTYPE'
!        RUC37_struc(n)%LDT_ncvar_albedo_monthly = 'RUC37_ALBEDO_MONTHLY'
!        RUC37_struc(n)%LDT_ncvar_shdfac_monthly = 'RUC37_SHDFAC_MONTHLY'
!        RUC37_struc(n)%LDT_ncvar_lai_monthly = 'RUC37_LAI_MONTHLY'
!        RUC37_struc(n)%LDT_ncvar_tbot        = 'RUC37_TBOT'
!        RUC37_struc(n)%LDT_ncvar_snoalb      = 'RUC37_SNOALB'

        RUC37_struc(n)%LDT_ncvar_vegetype    = 'LANDCOVER'
        RUC37_struc(n)%LDT_ncvar_soiltype    = 'TEXTURE'
        RUC37_struc(n)%LDT_ncvar_albedo_monthly = 'ALBEDO'
        RUC37_struc(n)%LDT_ncvar_shdfac_monthly = 'GREENNESS'
        RUC37_struc(n)%LDT_ncvar_lai_monthly = 'LAI'
        RUC37_struc(n)%LDT_ncvar_tbot        = 'TBOT'
        RUC37_struc(n)%LDT_ncvar_snoalb      = 'MXSNALBEDO'

    enddo

    ! set default restart format to netcdf
    do n=1,LIS_rc%nnest
        RUC37_struc(n)%rformat = "netcdf"
    enddo
    ! restart run, read restart file
    if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, "RUC37 restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, "RUC37 restart file: not defined")
        enddo
        
        Call ESMF_ConfigFindLabel(LIS_config, "RUC37 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, "RUC37 restart file format: not defined")
        enddo
    ! cold start run, read initial state variables
    else 
        ! surface emissivity (0.0 - 1.0).
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial emiss:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_emiss, rc=rc)
            call LIS_verify(rc, "RUC37 initial emiss: not defined")
        enddo

        ! exchange coefficient for head and moisture (m s-1).
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial ch:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_ch, rc=rc)
            call LIS_verify(rc, "RUC37 initial ch: not defined")
        enddo

        ! exchange coefficient for momentum (m s-1).
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial cm:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_cm, rc=rc)
            call LIS_verify(rc, "RUC37 initial cm: not defined")
        enddo

        ! water equivalent of accumulated snow depth (m).
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial sneqv:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_sneqv, rc=rc)
            call LIS_verify(rc, "RUC37 initial sneqv: not defined")
        enddo

        ! physical snow depth (m).
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial snowh:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_snowh, rc=rc)
            call LIS_verify(rc, "RUC37 initial snowh: not defined")
        enddo

        ! fractional snow cover ( fraction [0.0-1.0] )
        !call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial snowc:", rc = rc)
        do n=1,LIS_rc%nnest
            ! call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_snowc, rc=rc)
            ! call LIS_verify(rc, "RUC37 initial snowc: not defined")
            RUC37_struc(n)%init_snowc = min(1.0, RUC37_struc(n)%init_sneqv/0.016) 
        enddo
        
        ! canopy moisture content (kg m-2)
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial canwat:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_canwat, rc=rc)
            call LIS_verify(rc, "RUC37 initial canwat: not defined")
        enddo

        ! surface albedo including possible snow-cover effect.  this is set in lsmruc,
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial alb:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_alb, rc=rc)
            call LIS_verify(rc, "RUC37 initial alb: not defined")
        enddo

        ! total soil moisture content (m3 m-3)
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial smc:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, RUC37_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_smc(i), rc=rc)
            end do
            call LIS_verify(rc, "RUC37 initial smc: not defined")
        enddo

        ! liquid soil moisture content (m3 m-3)
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial sho:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, RUC37_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_sho(i), rc=rc)
            end do
            call LIS_verify(rc, "RUC37 initial sho: not defined")
        enddo

        ! soil temperature (k)
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial stc:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, RUC37_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_stc(i), rc=rc)
            end do
            call LIS_verify(rc, "RUC37 initial stc: not defined")
        enddo

        ! soil ice content (m3 m-3)
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial smfr:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, RUC37_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_smfr(i), rc=rc)
            end do
            call LIS_verify(rc, "RUC37 initial smfr: not defined")
        enddo

        ! frozen soil flag 
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial keepfr:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, RUC37_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_keepfr(i), rc=rc)
            end do
            call LIS_verify(rc, "RUC37 initial keepfr: not defined")
        enddo

        ! skin temperature (k)
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial tskin:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_tskin, rc=rc)
            call LIS_verify(rc, "RUC37 initial tskin: not defined")
        enddo

        ! mixing ratio at the surface ( kg kg{-1} )
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial qvg:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_qvg, rc=rc)
            call LIS_verify(rc, "RUC37 initial qvg: not defined")
        enddo

        ! sprcific humidity at the surface ( kg kg{-1} )
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial qsfc:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_qsfc, rc=rc)
            call LIS_verify(rc, "RUC37 initial qsfc: not defined")
        enddo

        ! cloud water mixing ratio at the surface ( kg kg{-1} )
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial qcg:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_qcg, rc=rc)
            call LIS_verify(rc, "RUC37 initial qcg: not defined")
        enddo

        ! surface water vapor mixing ratio at satration (kg kg-1)
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial qsg:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_qsg, rc=rc)
            call LIS_verify(rc, "RUC37 initial qsg: not defined")
        enddo

        ! snow temperature at 7.5 cm depth (k)
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial snt75cm:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_snt75cm, rc=rc)
            call LIS_verify(rc, "RUC37 initial snt75cm: not defined")
        enddo

        ! average snow temperature in k
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial tsnav:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_tsnav, rc=rc)
            call LIS_verify(rc, "RUC37 initial tsnav: not defined")
        enddo

        ! total soil column moisture content, frozen and unfrozen ( m )
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial soilm:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_soilm, rc=rc)
            call LIS_verify(rc, "RUC37 initial soilm: not defined")
        enddo

        ! available soil moisture in the root zone ( fraction [smcwlt-smcmax]
        call ESMF_ConfigFindLabel(LIS_config, "RUC37 initial smroot:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, RUC37_struc(n)%init_smroot, rc=rc)
            call LIS_verify(rc, "RUC37 initial smroot: not defined")
        enddo

    end if
     
    write(LIS_logunit, *) "Finish reading LIS configuration file for RUC37 model"
     
end subroutine RUC37_readcrd
