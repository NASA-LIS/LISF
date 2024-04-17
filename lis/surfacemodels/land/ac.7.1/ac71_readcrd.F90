!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: Ac71_readcrd
! \label{Ac71_readcrd}
!
! !REVISION HISTORY:
!   18 JAN 2024, Louise Busschaert; initial implementation for AC71
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
    character*32 :: soil_scheme_name

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
    

    ! First day of year of cropping period
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 starting day of crop period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%Crop_AnnualStartDay, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 starting day of crop period: not defined")
    enddo

    ! Last day of year of cropping period
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 ending day of crop period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%Crop_AnnualEndDay, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 ending day of crop period: not defined")
    enddo

    ! First month of year of cropping period
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 starting month of crop period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%Crop_AnnualStartMonth, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 starting month of crop period: not defined")
    enddo

    ! Last month of year of cropping period
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 ending month of crop period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%Crop_AnnualEndMonth, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 ending month of crop period: not defined")
    enddo

    ! number of soil compartments
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 max no of compartments:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%max_No_compartments, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 max no of compartments: not defined")
    enddo

    ! number of soil layers
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 number of soil layers:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%NrSoilLayers, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 number of soil layers: not defined")
    enddo
 
    ! allocate memory for sldpth using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(AC71_struc(n)%Thickness(AC71_struc(n)%NrSoilLayers))
        allocate(AC71_struc(n)%init_smc(AC71_struc(n)%NrSoilLayers))
    enddo

    ! PathNameSimul
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 input path:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
            AC71_struc(n)%PathNameSimul, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 input path: not defined")
    enddo

    ! CO2_Filename
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 CO2_Filename:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
            AC71_struc(n)%CO2_Filename, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 CO2_Filename: not defined")
    enddo
 
    ! Management_Filename
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 Management_Filename:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
            AC71_struc(n)%Management_Filename, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 Management_Filename: not defined")
    enddo
 
    ! Irrigation_Filename
    call ESMF_ConfigFindLabel(LIS_config, "Irrigation_Filename:", rc = rc)
    do n=1, LIS_rc%nnest
        if (rc == 0) then
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

 
    ! AquaCrop model soil parameter table
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 soil parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%soil_tbl_name, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.1 soil parameter table: not defined")
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
 
    do n=1, LIS_rc%nnest
      AC71_struc(n)%dt = AC71_struc(n)%ts
    enddo 

    ! reference height of forcings
    call ESMF_ConfigFindLabel(LIS_config, &
         "AquaCrop.7.1 reference height of forcings:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%refz_forc, rc=rc)
        call LIS_verify(rc, &
             "AquaCrop.7.1 reference height of forcings: "//&
             "not defined")
    enddo
 
    ! thickness of soil layers
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.1 soil layer thickness:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, AC71_struc(n)%NrSoilLayers
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%Thickness(i), rc=rc)
        call LIS_verify(rc, 'AquaCrop.7.1 soil layer thickness: not defined')
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
            call ESMF_ConfigGetAttribute(LIS_config, AC71_struc(n)%init_smc, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.1 initial liquid soil moistures: not defined")
        enddo
    endif

    deallocate(nids)

    write(LIS_logunit, *) "Finish reading LIS configuration file for AquaCrop.7.1 model"
end subroutine AC71_readcrd
