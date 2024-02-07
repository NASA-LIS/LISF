!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: Ac71_readcrd
! \label{Ac71_readcrd}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   9/4/14 : Shugong Wang, initial implementation for LIS 7 and Ac71
!   18 JAN 2024, Louise Busschaert; initial implementation for LIS 7 and AC71
!
! !INTERFACE:
subroutine Ac71_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit , LIS_verify, LIS_endrun
    use Ac71_lsmMod, only       : AC71_struc
    use netcdf 
!
! !DESCRIPTION:
!
!  This routine reads the options specific to Ac71 model from
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

    write(LIS_logunit, *) "Start reading LIS configuration file for AquaCrop.7.1 model"
    
    ! open NetCDF parameter file for reading global attributes 
    do n=1,LIS_rc%nnest
      ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE,ncid=nids(n))
      call LIS_verify(ios,'Error in nf90_open in '//trim(LIS_rc%paramfile(n))//' in Ac71_readcrd')
    enddo 
 
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "AquaCrop.7.1 model timestep: not defined")
        call LIS_parseTimeString(time, AC71_struc(n)%ts)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc,"AquaCrop.7.1 restart output interval: not defined")
        call LIS_parseTimeString(time, AC71_struc(n)%rstInterval)
    enddo
    
    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number of soil layers
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 number of soil layers:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%nsoil, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 number of soil layers: not defined")
    enddo
 
    ! MB: AC71

    ! First day of year of cropping period
    call ESMF_ConfigFindLabel(LIS_config, "Starting day of crop period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%Crop_AnnualStartDay, rc=rc)
        call LIS_verify(rc, "Starting day of crop period: not defined")
    enddo

    ! Last day of year of cropping period
    call ESMF_ConfigFindLabel(LIS_config, "Ending day of crop period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%Crop_AnnualEndDay, rc=rc)
        call LIS_verify(rc, "Ending day of crop period: not defined")
    enddo

    ! First month of year of cropping period
    call ESMF_ConfigFindLabel(LIS_config, "Starting month of crop period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%Crop_AnnualStartMonth, rc=rc)
        call LIS_verify(rc, "Starting month of crop period: not defined")
    enddo

    ! Last month of year of cropping period
    call ESMF_ConfigFindLabel(LIS_config, "Ending month of crop period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%Crop_AnnualEndMonth, rc=rc)
        call LIS_verify(rc, "Ending month of crop period: not defined")
    enddo

    ! number of soil layers
    call ESMF_ConfigFindLabel(LIS_config, "max_No_compartments:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%max_No_compartments, rc=rc)
        call LIS_verify(rc, "max_No_compartments: not defined")
    enddo

    ! number of soil layers
    call ESMF_ConfigFindLabel(LIS_config, "NrSoilLayers:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%NrSoilLayers, rc=rc)
        call LIS_verify(rc, "NrSoilLayers: not defined")
    enddo
 
    ! allocate memory for sldpth using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(AC71_struc(n)%Thickness(AC71_struc(n)%NrSoilLayers))
        allocate(AC71_struc(n)%sldpth(AC71_struc(n)%nsoil))
        allocate(AC71_struc(n)%init_stc( AC71_struc(n)%nsoil))
        allocate(AC71_struc(n)%init_sh2o(AC71_struc(n)%nsoil))
        allocate(AC71_struc(n)%init_smc(AC71_struc(n)%nsoil))
    enddo
 
    ! maximum number of snow layers
    do n=1, LIS_rc%nnest
       AC71_struc(n)%nsnow  = 3
    enddo
 
    ! MB: AC71
    ! PathNameOutp
    !call ESMF_ConfigFindLabel(LIS_config, "PathNameOutp:", rc = rc)
    !do n=1, LIS_rc%nnest
    !    call ESMF_ConfigGetAttribute(LIS_config, &
    !        AC71_struc(n)%PathNameOutp, rc=rc)
    !    call LIS_verify(rc, "PathNameOutp: not defined")
    !enddo
 
    ! PathNameSimul
    call ESMF_ConfigFindLabel(LIS_config, "AC_INPUT_PATH:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
            AC71_struc(n)%PathNameSimul, rc=rc)
        call LIS_verify(rc, "PathNameSimul: not defined")
    enddo
 
    ! PathNameList
    !call ESMF_ConfigFindLabel(LIS_config, "PathNameList:", rc = rc)
    !do n=1, LIS_rc%nnest
    !    call ESMF_ConfigGetAttribute(LIS_config, &
    !        AC71_struc(n)%PathNameList, rc=rc)
    !    call LIS_verify(rc, "PathNameList: not defined")
    !enddo
 
    ! PathNameParam
    !call ESMF_ConfigFindLabel(LIS_config, "PathNameParam:", rc = rc)
    !do n=1, LIS_rc%nnest
    !    call ESMF_ConfigGetAttribute(LIS_config, &
    !        AC71_struc(n)%PathNameParam, rc=rc)
    !    call LIS_verify(rc, "PathNameParam: not defined")
    !enddo

    ! CO2_Filename
    call ESMF_ConfigFindLabel(LIS_config, "CO2_Filename:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
            AC71_struc(n)%CO2_Filename, rc=rc)
        call LIS_verify(rc, "CO2_Filename: not defined")
    enddo
 
    ! Crop_Filename
    call ESMF_ConfigFindLabel(LIS_config, "Crop_Filename:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
            AC71_struc(n)%Crop_Filename, rc=rc)
        call LIS_verify(rc, "Crop_Filename: not defined")
    enddo
 
    ! Management_Filename
    call ESMF_ConfigFindLabel(LIS_config, "Management_Filename:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
            AC71_struc(n)%Management_Filename, rc=rc)
        call LIS_verify(rc, "Management_Filename: not defined")
    enddo
 
    ! Irrigation_Filename
    call ESMF_ConfigFindLabel(LIS_config, "Irrigation_Filename:", rc = rc)
    do n=1, LIS_rc%nnest
        if ( rc == 0) then
            call ESMF_ConfigGetAttribute(LIS_config, &
                 AC71_struc(n)%Irrigation_Filename, rc=rc)
             ! change lis none to AquaCrop (None)
             if ((AC71_struc(n)%Irrigation_Filename .eq. 'none') .or. &
                 (AC71_struc(n)%Irrigation_Filename .eq. 'None')) then
                 AC71_struc(n)%Irrigation_Filename = '(None)'
             endif 
        else
            print*,'Irrigation_Filename: not defined --> set to (None)'
            AC71_struc(n)%Irrigation_Filename = '(None)'
        endif
    enddo

 
    !! NumberSimulationRuns
    !call ESMF_ConfigFindLabel(LIS_config, "NumberSimulationRuns:", rc = rc)
    !do n=1, LIS_rc%nnest
    !    call ESMF_ConfigGetAttribute(LIS_config, &
    !        AC71_struc(n)%NumberSimulationRuns, rc=rc)
    !    call LIS_verify(rc, "NumberSimulationRuns: not defined")
    !enddo

    ! MB: AC71

 
    ! Noah model landuse parameter table
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 landuse parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%landuse_tbl_name, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 landuse parameter table: not defined")
    enddo
 
    ! Noah model soil parameter table
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 soil parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%soil_tbl_name, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 soil parameter table: not defined")
    enddo
 
    ! Noah model general parameter table
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 general parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%gen_tbl_name, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 general parameter table: not defined")
    enddo
 
    ! Ac parameter table
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 MP parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%ac_tbl_name, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 MP parameter table: not defined")
    enddo
 
    ! landuse classification scheme
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'LANDCOVER_SCHEME', landuse_scheme_name)
        call LIS_verify(ios, 'Error in nf90_get_att: LANDCOVER_SCHEME')
        if (trim(landuse_scheme_name) .eq. "USGS") then
          AC71_struc(n)%landuse_scheme_name = "USGS"
        elseif (trim(landuse_scheme_name) .eq. "IGBPNCEP") then
          AC71_struc(n)%landuse_scheme_name = "MODIFIED_IGBP_MODIS_NOAH"
        elseif (trim(landuse_scheme_name) .eq. "UMD") then
          AC71_struc(n)%landuse_scheme_name = "UMD"
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
          AC71_struc(n)%soil_scheme_name = "STAS"
        else
          write(LIS_logunit, *) "Fatal error: currently, only STATSGO soil scheme is supported by Noah-MP!"
          call LIS_endrun()
        endif 
    enddo
 
    ! vegetation model
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 vegetation model option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%dveg_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 vegetation model option: not defined")
    enddo
 
    ! canopy stomatal resistance
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 canopy stomatal resistance option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%crs_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 canopy stomatal resistance option: not defined")
    enddo
 
    ! soil moisture factor for stomatal resistance
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 soil moisture factor for stomatal resistance option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%btr_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 soil moisture factor for stomatal resistance option: not defined")
    enddo
 
    ! runoff and groundwater
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 runoff and groundwater option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%run_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 runoff and groundwater option: not defined")
    enddo
 
    ! surface layer drag coefficients (CH & CM)
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 surface layer drag coefficient option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%sfc_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 surface layer drag coefficient option: not defined")
    enddo
 
    ! supercooled liquid water
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 supercooled liquid water option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%frz_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 supercooled liquid water option: not defined")
    enddo
 
    ! frozen soil permeability
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 frozen soil permeability option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%inf_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 frozen soil permeability option: not defined")
    enddo
 
    ! radiation transfer
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 radiation transfer option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%rad_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 radiation transfer option: not defined")
    enddo
 
    ! snow surface albedo
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 snow surface albedo option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%alb_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 snow surface albedo option: not defined")
    enddo
 
    ! rainfall & snowfall
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 rainfall and snowfall option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%snf_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 rainfall and snowfall option: not defined")
    enddo
 
    ! lower boundary of soil temperature
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 lower boundary of soil temperature option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%tbot_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 lower boundary of soil temperature option: not defined")
    enddo
 
    ! snow/soil temperature time scheme
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 snow and soil temperature time scheme:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%stc_opt, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 snow and soil temperature time scheme: not defined")
    enddo
 
    ! the number of total soil types in parameter table
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'NUMBER_SOILTYPES', AC71_struc(n)%nslcats)
        call LIS_verify(ios, 'Error in nf90_get_att: NUMBER_SOILTYPES')
    enddo
 
    ! the number of total land cover types in parameter table
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'NUMBER_LANDCATS', AC71_struc(n)%nlucats)
        call LIS_verify(ios, 'Error in nf90_get_att: NUMBER_LANDCATS')
    enddo
 
    do n=1, LIS_rc%nnest
      AC71_struc(n)%dt = AC71_struc(n)%ts
    enddo 

    ! thickness of soil layers
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 soil layer thickness:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, AC71_struc(n)%nsoil
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%sldpth(i), rc=rc)
            call LIS_verify(rc, 'AquaCrop.7.1 soil layer thickness: not defined')
        enddo
    enddo
 
    ! MB: AC71
    ! thickness of soil layers
    call ESMF_ConfigFindLabel(LIS_config, "Thickness:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, AC71_struc(n)%NrSoilLayers
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%Thickness(i), rc=rc)
            call LIS_verify(rc, 'Thickness: not defined')
        enddo
    enddo
 
    ! urban land cover type index
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'URBANCLASS', AC71_struc(n)%urban_vegetype)
        call LIS_verify(ios, 'Error in nf90_get_att: URBANCLASS')
    enddo
 
    ! ice flag: 0 = no ice, 1 = ice
    do n=1, LIS_rc%nnest
        AC71_struc(n)%ice_flag = 0
    enddo
 
    ! surface type 1=soil, 2=lake
    do n=1, LIS_rc%nnest
        AC71_struc(n)%st_flag = 1 
    enddo
 
    ! soil color type
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 soil color index:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%sc_idx, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 soil color index: not defined")
    enddo
 
    ! option of Chen adjustment of Czil
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 CZIL option (iz0tlnd):", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%iz0tlnd, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 CZIL option (iz0tlnd): not defined")
    enddo
 
    do n=1,LIS_rc%nnest
      ios = nf90_close(nids(n))
      call LIS_verify(ios,'Error in nf90_close in '//trim(LIS_rc%paramfile(n))//' in Ac71_readcrd')
    enddo 

    ! The following lines hard code the LDT NetCDF variable names. 
    do n=1, LIS_rc%nnest
        AC71_struc(n)%LDT_ncvar_shdfac_monthly = 'GREENNESS'  !'AC71_SHDFAC_MONTHLY'
        ! AC71_struc(n)%LDT_ncvar_vegetype = ' ! Edit here if hard code name
        AC71_struc(n)%LDT_ncvar_soiltype = 'AC71_SOILTYPE'
        AC71_struc(n)%LDT_ncvar_slopetype = 'SLOPETYPE'       !'AC71_SLOPETYPE'
        AC71_struc(n)%LDT_ncvar_smceq    = 'AC71_SMCEQ'
        AC71_struc(n)%LDT_ncvar_tbot     = 'TBOT'             !'AC71_TBOT'
        AC71_struc(n)%LDT_ncvar_pblh     = 'NOAHMP36_PBLH'
    enddo

    ! set default restart format to netcdf
    do n=1,LIS_rc%nnest
        AC71_struc(n)%rformat = "netcdf"
    enddo
    ! restart run, read restart file
    if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 restart file: not defined")
        enddo
        
        Call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 restart file format: not defined")
        enddo
    ! cold start run, read initial state variables
    else 
        ! snow albedo at last time step
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial value of snow albedo at the last timestep:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_albold, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial value of snow albedo at the last timestep: not defined")
        enddo

        ! snow mass at the last time step
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial value of snow mass at the last timestep:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_sneqvo, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial value of snow mass at the last timestep: not defined")
        enddo

        ! soil temperature
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial soil temperatures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, AC71_struc(n)%nsoil 
                call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_stc(i), rc=rc)
            end do
            call LIS_verify(rc, "AquaCrop.7.1 initial soil temperatures: not defined")
        enddo

        ! volumetric liquid soil moisture
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial liquid soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, AC71_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_sh2o(i), rc=rc)
            end do
            call LIS_verify(rc, "AquaCrop.7.1 initial liquid soil moistures: not defined")
        enddo

        ! volumetric soil moisture, ice + liquid
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial total soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, AC71_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_smc(i), rc=rc)
            end do
            call LIS_verify(rc, "AquaCrop.7.1 initial total soil moistures: not defined")
        enddo

        ! canopy air temperature
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial canopy air temperature:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_tah, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial canopy air temperature: not defined")
        enddo

        ! canopy air vapor pressure
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial canopy air vapor pressure:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_eah, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial canopy air vapor pressure: not defined")
        enddo

        ! wetted or snowed fraction of canopy
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial wetted or snowed fraction of canopy:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_fwet, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial wetted or snowed fraction of canopy: not defined")
        enddo

        ! intercepted liquid water
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial intercepted liquid water:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_canliq, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial intercepted liquid water: not defined")
        enddo

        ! intercepted ice mass
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial intercepted ice mass:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_canice, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial intercepted ice mass: not defined")
        enddo

        ! vegetation temperature
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial vegetation temperature:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_tv, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial vegetation temperature: not defined")
        enddo

        ! ground temperature (skin temperature)
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial ground temperature:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_tg, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial ground temperature: not defined")
        enddo

        ! snowfall on the ground
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial snowfall on the ground:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_qsnow, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial snowfall on the ground: not defined")
        enddo

        ! snow height
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial snow height:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_snowh, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial snow height: not defined")
        enddo

        ! snow water equivalent
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial snow water equivalent:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_sneqv, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial snow water equivalent: not defined")
        enddo


        ! depth to water table
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial depth to water table:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_zwt, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial depth to water table: not defined")
        enddo

        ! water storage in aquifer
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial water storage in aquifer:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_wa, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial water storage in aquifer: not defined")
        enddo

        ! water in aquifer and saturated soil
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial water in aquifer and saturated soil:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_wt, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial water in aquifer and saturated soil: not defined")
        enddo

        ! lake water storage
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial lake water storage:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_wslake, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial lake water storage: not defined")
        enddo

        ! leaf mass
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial leaf mass:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_lfmass, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial leaf mass: not defined")
        enddo

        ! mass of fine roots
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial mass of fine roots:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_rtmass, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial mass of fine roots: not defined")
        enddo

        ! stem mass
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial stem mass:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_stmass, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial stem mass: not defined")
        enddo

        ! mass of wood including woody roots
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial mass of wood including woody roots:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_wood, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial mass of wood including woody roots: not defined")
        enddo

        ! stable carbon in deep soil
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial stable carbon in deep soil:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_stblcp, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial stable carbon in deep soil: not defined")
        enddo

        ! short-lived carbon in shallow soil
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial short-lived carbon in shallow soil:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_fastcp, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial short-lived carbon in shallow soil: not defined")
        enddo

        ! leaf area index
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial LAI:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_lai, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial LAI: not defined")
        enddo

        ! stem area index
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial SAI:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_sai, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial SAI: not defined")
        enddo

        ! momentum drag coefficient
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial momentum drag coefficient:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_cm, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial momentum drag coefficient: not defined")
        enddo

        ! sensible heat exchange coefficient
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial sensible heat exchange coefficient:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_ch, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial sensible heat exchange coefficient: not defined")
        enddo

        ! snow aging term
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial snow aging term:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_tauss, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial snow aging term: not defined")
        enddo

        ! soil water content between bottom of the soil and water table
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial soil water content between bottom of the soil and water table:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_smcwtd, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial soil water content between bottom of the soil and water table: not defined")
        enddo

        ! recharge to or from the water table when deep
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial recharge to or from the water table when deep:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_deeprech, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial recharge to or from the water table when deep: not defined")
        enddo

        ! recharge to or from the water table when shallow
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial recharge to or from the water table when shallow:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_rech, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial recharge to or from the water table when shallow: not defined")
        enddo
        
        ! air temperature and humidity reference height
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial reference height of temperature and humidity:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_zlvl, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial reference height of temperature and humidity: not defined")
        enddo

    end if
     
    deallocate(nids)

    write(LIS_logunit, *) "Finish reading LIS configuration file for AquaCrop.7.1 model"
end subroutine AC71_readcrd
