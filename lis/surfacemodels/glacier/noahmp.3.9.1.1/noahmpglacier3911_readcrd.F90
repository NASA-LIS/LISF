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
! !ROUTINE: noahmpglacier3911_readcrd
! \label{noahmpglacier3911_readcrd}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine noahmpglacier3911_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit , LIS_verify, LIS_endrun
    use noahmpglacier3911_Mod, only       : noahmpgl3911_struc
    use netcdf 
!
! !DESCRIPTION:
!
!  This routine reads the options specific to NoahMP36 model from
!  the LIS configuration file.
!
!EOP
    implicit none

    integer      :: rc 
    integer      :: n, i
    character*10 :: time 
    character*6  :: str_i
    integer :: ios
    character*32 :: soil_scheme_name, landuse_scheme_name

    write(LIS_logunit, *) "[INFO] Start reading LIS configuration file for Noah-MP glacier model version 3.9.1.1"
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 model timestep: not defined")
        call LIS_parseTimeString(time, Noahmpgl3911_struc(n)%ts)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc,"Noah-MP glacier 3.9.1.1 restart output interval: not defined")
        call LIS_parseTimeString(time, Noahmpgl3911_struc(n)%rstInterval)
    enddo
    
! snow surface albedo
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 snow surface albedo option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, noahmpgl3911_struc(n)%alb_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 snow surface albedo option: not defined")
    enddo

    ! rainfall & snowfall
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 rainfall and snowfall option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, noahmpgl3911_struc(n)%snf_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 rainfall and snowfall option: not defined")
    enddo
 
    ! lower boundary of soil temperature
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 lower boundary of soil temperature option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, noahmpgl3911_struc(n)%tbot_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 lower boundary of soil temperature option: not defined")
    enddo
 
    ! snow/soil temperature time scheme
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 snow and soil temperature time scheme:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, noahmpgl3911_struc(n)%stc_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 snow and soil temperature time scheme: not defined")
    enddo

    ! glacier treatment option
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 glacier treatment scheme:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, noahmpgl3911_struc(n)%gla_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 glacier treatment scheme: not defined")
    enddo
    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number of soil layersn
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 number of soil layers:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%nsoil, rc=rc)
        call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 number of soil layers: not defined")
    enddo
 
    ! allocate memory for sldpth using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(Noahmpgl3911_struc(n)%sldpth(Noahmpgl3911_struc(n)%nsoil))
        allocate(Noahmpgl3911_struc(n)%init_stc( Noahmpgl3911_struc(n)%nsoil))
        allocate(Noahmpgl3911_struc(n)%init_sh2o(Noahmpgl3911_struc(n)%nsoil))
        allocate(Noahmpgl3911_struc(n)%init_smc(Noahmpgl3911_struc(n)%nsoil))
    enddo
 
    ! maximum number of snow layers
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 number of snow layers:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%nsnow, rc=rc)
        call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 number of snow layers: not defined")
    enddo
 
    ! air temperature and humidity reference height
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial reference height of temperature and humidity:", rc = rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_zlvl, rc=rc)
       call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial reference height of temperature and humidity: not defined")
    enddo

    ! snowfall on the ground
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial snowfall on the ground:", rc = rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_qsnow, rc=rc)
       call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial snowfall on the ground: not defined")
    enddo

    ! snow mass at the last time step
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial value of snow mass at the last timestep:", rc = rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_sneqvo, rc=rc)
       call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial value of snow mass at the last timestep: not defined")
    enddo
 ! snow height
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial snow height:", rc = rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_snowh, rc=rc)
       call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial snow height: not defined")
    enddo

        ! snow water equivalent
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial snow water equivalent:", rc = rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_sneqv, rc=rc)
       call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial snow water equivalent: not defined")
    enddo

        ! ground temperature (skin temperature)
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial ground temperature:", rc = rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_tg, rc=rc)
       call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial ground temperature: not defined")
    enddo
    
  ! The following lines hard code the LDT NetCDF variable names. 
    do n=1, LIS_rc%nnest
        noahmpgl3911_struc(n)%LDT_ncvar_shdfac_monthly = 'GREENNESS'  !'NOAHMP36_SHDFAC_MONTHLY'
        ! Noahmpgl3911_struc(n)%LDT_ncvar_vegetype = ' ! Edit here if hard code name
        noahmpgl3911_struc(n)%LDT_ncvar_soiltype = 'NOAHMP36_SOILTYPE'
        noahmpgl3911_struc(n)%LDT_ncvar_slopetype = 'SLOPETYPE'       !'NOAHMP36_SLOPETYPE'
        noahmpgl3911_struc(n)%LDT_ncvar_smceq    = 'NOAHMP36_SMCEQ'
        noahmpgl3911_struc(n)%LDT_ncvar_tbot     = 'TBOT'             !'NOAHMP36_TBOT'
        noahmpgl3911_struc(n)%LDT_ncvar_pblh     = 'NOAHMP36_PBLH'
    enddo

    ! set default restart format to netcdf
    do n=1,LIS_rc%nnest
        Noahmpgl3911_struc(n)%rformat = "netcdf"
     enddo
     ! thickness of soil layers
     call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 soil layer thickness:", rc = rc)
     do n=1, LIS_rc%nnest
        do i = 1, Noahmpgl3911_struc(n)%nsoil
           call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%sldpth(i), rc=rc)
           call LIS_verify(rc, 'Noah-MP glacier 3.9.1.1 soil layer thickness: not defined')
        enddo
     enddo
    ! restart run, read restart file
     if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 restart file: not defined")
        enddo
        
        Call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 restart file format: not defined")
        enddo
    ! cold start run, read initial state variables
    else 
        ! snow albedo at last time step
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial value of snow albedo at the last timestep:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_albold, rc=rc)
            call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial value of snow albedo at the last timestep: not defined")
        enddo

        ! momentum drag coefficient
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial momentum drag coefficient:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_cm, rc=rc)
            call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial momentum drag coefficient: not defined")
        enddo

        ! sensible heat exchange coefficient
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial sensible heat exchange coefficient:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_ch, rc=rc)
            call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial sensible heat exchange coefficient: not defined")
        enddo

        ! snow aging term
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial snow aging term:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_tauss, rc=rc)
            call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial snow aging term: not defined")
        enddo

 ! soil temperature
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial soil temperatures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, Noahmpgl3911_struc(n)%nsoil 
                call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_stc(i), rc=rc)
            end do
            call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial soil temperatures: not defined")
        enddo

        ! volumetric liquid soil moisture
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial liquid soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, Noahmpgl3911_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_sh2o(i), rc=rc)
            end do
            call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial liquid soil moistures: not defined")
        enddo

        ! volumetric soil moisture, ice + liquid
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP glacier 3.9.1.1 initial total soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, Noahmpgl3911_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, Noahmpgl3911_struc(n)%init_smc(i), rc=rc)
            end do
            call LIS_verify(rc, "Noah-MP glacier 3.9.1.1 initial total soil moistures: not defined")
        enddo


     end if
    write(LIS_logunit, *) "Finish reading LIS configuration file for Noah-MP glacier 3.9.1.1 model"

  end subroutine Noahmpglacier3911_readcrd
