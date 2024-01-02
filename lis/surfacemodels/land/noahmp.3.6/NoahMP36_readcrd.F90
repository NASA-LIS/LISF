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
! !ROUTINE: NoahMP36_readcrd
! \label{NoahMP36_readcrd}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   9/4/14 : Shugong Wang, initial implementation for LIS 7 and NoahMP36
!
! !INTERFACE:
subroutine NoahMP36_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit , LIS_verify, LIS_endrun
    use NoahMP36_lsmMod, only       : NOAHMP36_struc
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
    integer, allocatable :: nids(:)
    character*32 :: soil_scheme_name, landuse_scheme_name

    allocate(nids(LIS_rc%nnest))

    write(LIS_logunit, *) "Start reading LIS configuration file for Noah-MP.3.6 model"
    
    ! open NetCDF parameter file for reading global attributes 
    do n=1,LIS_rc%nnest
      ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE,ncid=nids(n))
      call LIS_verify(ios,'Error in nf90_open in '//trim(LIS_rc%paramfile(n))//' in NoahMP36_readcrd')
    enddo 
 
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "Noah-MP.3.6 model timestep: not defined")
        call LIS_parseTimeString(time, NOAHMP36_struc(n)%ts)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc,"Noah-MP.3.6 restart output interval: not defined")
        call LIS_parseTimeString(time, NOAHMP36_struc(n)%rstInterval)
    enddo
    
    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number of soil layers
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 number of soil layers:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%nsoil, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 number of soil layers: not defined")
    enddo
 
    ! allocate memory for sldpth using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(NOAHMP36_struc(n)%sldpth(NOAHMP36_struc(n)%nsoil))
        allocate(NOAHMP36_struc(n)%init_stc( NOAHMP36_struc(n)%nsoil))
        allocate(NOAHMP36_struc(n)%init_sh2o(NOAHMP36_struc(n)%nsoil))
        allocate(NOAHMP36_struc(n)%init_smc(NOAHMP36_struc(n)%nsoil))
    enddo
 
    ! maximum number of snow layers
    do n=1, LIS_rc%nnest
       NOAHMP36_struc(n)%nsnow  = 3
    enddo
 
 
    ! Noah model landuse parameter table
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 landuse parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%landuse_tbl_name, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 landuse parameter table: not defined")
    enddo
 
    ! Noah model soil parameter table
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 soil parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%soil_tbl_name, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 soil parameter table: not defined")
    enddo
 
    ! Noah model general parameter table
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 general parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%gen_tbl_name, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 general parameter table: not defined")
    enddo
 
    ! NoahMP parameter table
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 MP parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%noahmp_tbl_name, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 MP parameter table: not defined")
    enddo
 
    ! landuse classification scheme
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'LANDCOVER_SCHEME', landuse_scheme_name)
        call LIS_verify(ios, 'Error in nf90_get_att: LANDCOVER_SCHEME')
        if (trim(landuse_scheme_name) .eq. "USGS") then
          NOAHMP36_struc(n)%landuse_scheme_name = "USGS"
        elseif (trim(landuse_scheme_name) .eq. "IGBPNCEP") then
          NOAHMP36_struc(n)%landuse_scheme_name = "MODIFIED_IGBP_MODIS_NOAH"
        elseif (trim(landuse_scheme_name) .eq. "UMD") then
          NOAHMP36_struc(n)%landuse_scheme_name = "UMD"
        else
          write(LIS_logunit, *) "Fatal error: currently, only USGS and IGBPNCEP is supported by Noah-MP!"
          call LIS_endrun()
        endif
    enddo
 
    ! soil classification scheme
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'SOILTEXT_SCHEME', soil_scheme_name)
        call LIS_verify(ios, 'Error in nf90_get_att: SOILTEXT_SCHEME')
        if (trim(soil_scheme_name) .eq. "STATSGO") then 
          NOAHMP36_struc(n)%soil_scheme_name = "STAS"
        else
          write(LIS_logunit, *) "Fatal error: currently, only STATSGO soil scheme is supported by Noah-MP!"
          call LIS_endrun()
        endif 
    enddo
 
    ! vegetation model
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 vegetation model option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%dveg_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 vegetation model option: not defined")
    enddo
 
    ! canopy stomatal resistance
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 canopy stomatal resistance option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%crs_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 canopy stomatal resistance option: not defined")
    enddo
 
    ! soil moisture factor for stomatal resistance
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 soil moisture factor for stomatal resistance option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%btr_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 soil moisture factor for stomatal resistance option: not defined")
    enddo
 
    ! runoff and groundwater
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 runoff and groundwater option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%run_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 runoff and groundwater option: not defined")
    enddo
 
    ! surface layer drag coefficients (CH & CM)
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 surface layer drag coefficient option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%sfc_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 surface layer drag coefficient option: not defined")
    enddo
 
    ! supercooled liquid water
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 supercooled liquid water option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%frz_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 supercooled liquid water option: not defined")
    enddo
 
    ! frozen soil permeability
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 frozen soil permeability option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%inf_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 frozen soil permeability option: not defined")
    enddo
 
    ! radiation transfer
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 radiation transfer option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%rad_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 radiation transfer option: not defined")
    enddo
 
    ! snow surface albedo
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 snow surface albedo option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%alb_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 snow surface albedo option: not defined")
    enddo
 
    ! rainfall & snowfall
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 rainfall and snowfall option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%snf_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 rainfall and snowfall option: not defined")
    enddo
 
    ! lower boundary of soil temperature
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 lower boundary of soil temperature option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%tbot_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 lower boundary of soil temperature option: not defined")
    enddo
 
    ! snow/soil temperature time scheme
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 snow and soil temperature time scheme:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%stc_opt, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 snow and soil temperature time scheme: not defined")
    enddo
 
    ! the number of total soil types in parameter table
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'NUMBER_SOILTYPES', NOAHMP36_struc(n)%nslcats)
        call LIS_verify(ios, 'Error in nf90_get_att: NUMBER_SOILTYPES')
    enddo
 
    ! the number of total land cover types in parameter table
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'NUMBER_LANDCATS', NOAHMP36_struc(n)%nlucats)
        call LIS_verify(ios, 'Error in nf90_get_att: NUMBER_LANDCATS')
    enddo
 
    ! the number of total slope category for Noah baseflow
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'NUMBER_SLOPETYPES', NOAHMP36_struc(n)%nslpcats)
        call LIS_verify(ios, 'Error in nf90_get_att: NUMBER_SLOPETYPES')
    enddo
 
    do n=1, LIS_rc%nnest
      NOAHMP36_struc(n)%dt = NOAHMP36_struc(n)%ts
    enddo 

    ! thickness of soil layers
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 soil layer thickness:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, NOAHMP36_struc(n)%nsoil
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%sldpth(i), rc=rc)
            call LIS_verify(rc, 'Noah-MP.3.6 soil layer thickness: not defined')
        enddo
    enddo
 
    ! urban land cover type index
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'URBANCLASS', NOAHMP36_struc(n)%urban_vegetype)
        call LIS_verify(ios, 'Error in nf90_get_att: URBANCLASS')
    enddo
 
    ! ice flag: 0 = no ice, 1 = ice
    do n=1, LIS_rc%nnest
        NOAHMP36_struc(n)%ice_flag = 0
    enddo
 
    ! surface type 1=soil, 2=lake
    do n=1, LIS_rc%nnest
        NOAHMP36_struc(n)%st_flag = 1 
    enddo
 
    ! soil color type
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 soil color index:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%sc_idx, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 soil color index: not defined")
    enddo
 
    ! option of Chen adjustment of Czil
    call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 CZIL option (iz0tlnd):", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%iz0tlnd, rc=rc)
        call LIS_verify(rc, "Noah-MP.3.6 CZIL option (iz0tlnd): not defined")
    enddo
 
    do n=1,LIS_rc%nnest
      ios = nf90_close(nids(n))
      call LIS_verify(ios,'Error in nf90_close in '//trim(LIS_rc%paramfile(n))//' in NoahMP36_readcrd')
    enddo 

    ! The following lines hard code the LDT NetCDF variable names. 
    do n=1, LIS_rc%nnest
        NOAHMP36_struc(n)%LDT_ncvar_shdfac_monthly = 'GREENNESS'  !'NOAHMP36_SHDFAC_MONTHLY'
        ! NOAHMP36_struc(n)%LDT_ncvar_vegetype = ' ! Edit here if hard code name
        NOAHMP36_struc(n)%LDT_ncvar_soiltype = 'NOAHMP36_SOILTYPE'
        NOAHMP36_struc(n)%LDT_ncvar_slopetype = 'SLOPETYPE'       !'NOAHMP36_SLOPETYPE'
        NOAHMP36_struc(n)%LDT_ncvar_smceq    = 'NOAHMP36_SMCEQ'
        NOAHMP36_struc(n)%LDT_ncvar_tbot     = 'TBOT'             !'NOAHMP36_TBOT'
        NOAHMP36_struc(n)%LDT_ncvar_pblh     = 'NOAHMP36_PBLH'
    enddo

    ! set default restart format to netcdf
    do n=1,LIS_rc%nnest
        NOAHMP36_struc(n)%rformat = "netcdf"
    enddo
    ! restart run, read restart file
    if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 restart file: not defined")
        enddo
        
        Call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 restart file format: not defined")
        enddo
    ! cold start run, read initial state variables
    else 
        ! snow albedo at last time step
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial value of snow albedo at the last timestep:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_albold, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial value of snow albedo at the last timestep: not defined")
        enddo

        ! snow mass at the last time step
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial value of snow mass at the last timestep:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_sneqvo, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial value of snow mass at the last timestep: not defined")
        enddo

        ! soil temperature
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial soil temperatures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, NOAHMP36_struc(n)%nsoil 
                call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_stc(i), rc=rc)
            end do
            call LIS_verify(rc, "Noah-MP.3.6 initial soil temperatures: not defined")
        enddo

        ! volumetric liquid soil moisture
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial liquid soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, NOAHMP36_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_sh2o(i), rc=rc)
            end do
            call LIS_verify(rc, "Noah-MP.3.6 initial liquid soil moistures: not defined")
        enddo

        ! volumetric soil moisture, ice + liquid
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial total soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, NOAHMP36_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_smc(i), rc=rc)
            end do
            call LIS_verify(rc, "Noah-MP.3.6 initial total soil moistures: not defined")
        enddo

        ! canopy air temperature
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial canopy air temperature:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_tah, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial canopy air temperature: not defined")
        enddo

        ! canopy air vapor pressure
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial canopy air vapor pressure:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_eah, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial canopy air vapor pressure: not defined")
        enddo

        ! wetted or snowed fraction of canopy
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial wetted or snowed fraction of canopy:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_fwet, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial wetted or snowed fraction of canopy: not defined")
        enddo

        ! intercepted liquid water
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial intercepted liquid water:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_canliq, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial intercepted liquid water: not defined")
        enddo

        ! intercepted ice mass
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial intercepted ice mass:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_canice, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial intercepted ice mass: not defined")
        enddo

        ! vegetation temperature
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial vegetation temperature:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_tv, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial vegetation temperature: not defined")
        enddo

        ! ground temperature (skin temperature)
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial ground temperature:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_tg, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial ground temperature: not defined")
        enddo

        ! snowfall on the ground
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial snowfall on the ground:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_qsnow, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial snowfall on the ground: not defined")
        enddo

        ! snow height
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial snow height:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_snowh, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial snow height: not defined")
        enddo

        ! snow water equivalent
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial snow water equivalent:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_sneqv, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial snow water equivalent: not defined")
        enddo


        ! depth to water table
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial depth to water table:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_zwt, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial depth to water table: not defined")
        enddo

        ! water storage in aquifer
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial water storage in aquifer:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_wa, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial water storage in aquifer: not defined")
        enddo

        ! water in aquifer and saturated soil
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial water in aquifer and saturated soil:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_wt, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial water in aquifer and saturated soil: not defined")
        enddo

        ! lake water storage
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial lake water storage:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_wslake, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial lake water storage: not defined")
        enddo

        ! leaf mass
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial leaf mass:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_lfmass, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial leaf mass: not defined")
        enddo

        ! mass of fine roots
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial mass of fine roots:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_rtmass, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial mass of fine roots: not defined")
        enddo

        ! stem mass
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial stem mass:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_stmass, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial stem mass: not defined")
        enddo

        ! mass of wood including woody roots
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial mass of wood including woody roots:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_wood, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial mass of wood including woody roots: not defined")
        enddo

        ! stable carbon in deep soil
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial stable carbon in deep soil:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_stblcp, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial stable carbon in deep soil: not defined")
        enddo

        ! short-lived carbon in shallow soil
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial short-lived carbon in shallow soil:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_fastcp, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial short-lived carbon in shallow soil: not defined")
        enddo

        ! leaf area index
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial LAI:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_lai, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial LAI: not defined")
        enddo

        ! stem area index
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial SAI:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_sai, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial SAI: not defined")
        enddo

        ! momentum drag coefficient
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial momentum drag coefficient:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_cm, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial momentum drag coefficient: not defined")
        enddo

        ! sensible heat exchange coefficient
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial sensible heat exchange coefficient:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_ch, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial sensible heat exchange coefficient: not defined")
        enddo

        ! snow aging term
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial snow aging term:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_tauss, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial snow aging term: not defined")
        enddo

        ! soil water content between bottom of the soil and water table
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial soil water content between bottom of the soil and water table:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_smcwtd, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial soil water content between bottom of the soil and water table: not defined")
        enddo

        ! recharge to or from the water table when deep
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial recharge to or from the water table when deep:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_deeprech, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial recharge to or from the water table when deep: not defined")
        enddo

        ! recharge to or from the water table when shallow
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial recharge to or from the water table when shallow:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_rech, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial recharge to or from the water table when shallow: not defined")
        enddo
        
        ! air temperature and humidity reference height
        call ESMF_ConfigFindLabel(LIS_config, "Noah-MP.3.6 initial reference height of temperature and humidity:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP36_struc(n)%init_zlvl, rc=rc)
            call LIS_verify(rc, "Noah-MP.3.6 initial reference height of temperature and humidity: not defined")
        enddo

    end if
     
    deallocate(nids)

    write(LIS_logunit, *) "Finish reading LIS configuration file for Noah-MP.3.6 model"
end subroutine NOAHMP36_readcrd
