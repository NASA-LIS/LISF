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
        allocate(AC71_struc(n)%init_smc(AC71_struc(n)%nsoil))
    enddo

    ! PathNameSimul
    call ESMF_ConfigFindLabel(LIS_config, "AC_INPUT_PATH:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
            AC71_struc(n)%PathNameSimul, rc=rc)
        call LIS_verify(rc, "PathNameSimul: not defined")
    enddo

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

 
    ! Noah model landuse parameter table
!    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 landuse parameter table:", rc = rc)
!    do n=1, LIS_rc%nnest
!        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%landuse_tbl_name, rc=rc)
!        call LIS_verify(rc, "AquaCrop.7.1 landuse parameter table: not defined")
!    enddo
 
    ! Noah model soil parameter table
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 soil parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%soil_tbl_name, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 soil parameter table: not defined")
    enddo
 
    ! Noah model general parameter table
!    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 general parameter table:", rc = rc)
!    do n=1, LIS_rc%nnest
!        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%gen_tbl_name, rc=rc)
!        call LIS_verify(rc, "AquaCrop.7.1 general parameter table: not defined")
!    enddo

 
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
        if (trim(soil_scheme_name) .eq. "STATSGO") then !LB: later change it to AC
          AC71_struc(n)%soil_scheme_name = "STAS"
        else
          write(LIS_logunit, *) "Fatal error: currently, only STATSGO soil scheme is supported by AquaCrop!"
          call LIS_endrun()
        endif 
    enddo
 
    ! the number of total land cover types in parameter table
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'NUMBER_LANDCATS', AC71_struc(n)%nlucats)
        call LIS_verify(ios, 'Error in nf90_get_att: NUMBER_LANDCATS')
    enddo
 
    do n=1, LIS_rc%nnest
      AC71_struc(n)%dt = AC71_struc(n)%ts
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

 
    do n=1,LIS_rc%nnest
      ios = nf90_close(nids(n))
      call LIS_verify(ios,'Error in nf90_close in '//trim(LIS_rc%paramfile(n))//' in Ac71_readcrd')
    enddo 

    ! The following lines hard code the LDT NetCDF variable names. 
    do n=1, LIS_rc%nnest
        AC71_struc(n)%LDT_ncvar_soiltype = 'AC71_SOILTYPE'
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
        ! volumetric liquid soil moisture
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 initial liquid soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            call LIS_verify(rc, "AquaCrop.7.1 initial liquid soil moistures: not defined")
        enddo
    endif

    deallocate(nids)

    write(LIS_logunit, *) "Finish reading LIS configuration file for AquaCrop.7.1 model"
end subroutine AC71_readcrd
