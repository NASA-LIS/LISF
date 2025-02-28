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
! !ROUTINE: AC72_readcrd
! \label{AC72_readcrd}
!
! !REVISION HISTORY:
!   04 NOV 2024, Louise Busschaert; initial implementation for AC72
!
! !INTERFACE:
subroutine AC72_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit , LIS_verify, LIS_endrun
    use AC72_lsmMod, only       : AC72_struc
    use netcdf 
!
! !DESCRIPTION:
!
!  This routine reads the options specific to AC72 model from
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
    character*32 :: soil_scheme_name,str

    allocate(nids(LIS_rc%nnest))

    write(LIS_logunit, *) "[INFO] Start reading LIS configuration file for AquaCrop.7.2 model"
    
    ! open NetCDF parameter file for reading global attributes 
    do n=1,LIS_rc%nnest
      ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE,ncid=nids(n))
      call LIS_verify(ios,'Error in nf90_open in '//trim(LIS_rc%paramfile(n))//' in AC72_readcrd')
    enddo 
 
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "AquaCrop.7.2 model timestep: not defined")
        call LIS_parseTimeString(time, AC72_struc(n)%ts)
        if (AC72_struc(n)%ts.ne.86400) then
          write(LIS_logunit, *) "Fatal error: AquaCrop.7.2 only runs with a daily time step ..."
          write(LIS_logunit, *) "Please select AquaCrop.7.2 model timestep: 1da"
          call LIS_endrun()
        endif
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc,"AquaCrop.7.2 restart output interval: not defined")
        call LIS_parseTimeString(time, AC72_struc(n)%rstInterval)
    enddo

    ! First day of year of simulation period
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 starting day of sim period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC72_struc(n)%Sim_AnnualStartDay, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.2 starting day of sim period: not defined")
    enddo

    ! First month of year of simulation period
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 starting month of sim period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC72_struc(n)%Sim_AnnualStartMonth, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.2 starting month of sim period: not defined")
    enddo

    ! First day of year of cropping period
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 starting day of crop period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC72_struc(n)%Crop_AnnualStartDay, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.2 starting day of crop period: not defined")
    enddo


    ! First month of year of cropping period
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 starting month of crop period:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC72_struc(n)%Crop_AnnualStartMonth, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.2 starting month of crop period: not defined")
    enddo

    ! PathNameSimul
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 input path:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
            AC72_struc(n)%PathNameSimul, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.2 input path: not defined")
    enddo

    ! CO2_Filename
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 CO2_Filename:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
            AC72_struc(n)%CO2_Filename, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.2 CO2_Filename: not defined")
    enddo
 
    ! Management_Filename
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 Management_Filename:", rc = rc)
    do n=1, LIS_rc%nnest
        if (rc == 0) then
            call ESMF_ConfigGetAttribute(LIS_config, &
                AC72_struc(n)%Management_Filename, rc=rc)
             ! change lis none to AquaCrop (None)
             if ((AC72_struc(n)%Management_Filename .eq. 'none') .or. &
                 (AC72_struc(n)%Management_Filename .eq. 'None')) then
                 AC72_struc(n)%Management_Filename = '(None)'
             endif 
        else
            write(LIS_logunit, *)'[INFO] AC72 Management_Filename: not defined, default management'
            AC72_struc(n)%Management_Filename = '(None)'
        endif
    enddo
 
    ! Irrigation_Filename
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 Irrigation_Filename:", rc = rc)
    do n=1, LIS_rc%nnest
        if (rc == 0) then
            call ESMF_ConfigGetAttribute(LIS_config, &
                 AC72_struc(n)%Irrigation_Filename, rc=rc)
             ! change lis none to AquaCrop (None)
             if ((AC72_struc(n)%Irrigation_Filename .eq. 'none') .or. &
                 (AC72_struc(n)%Irrigation_Filename .eq. 'None')) then
                 AC72_struc(n)%Irrigation_Filename = '(None)'
             endif 
        else
            write(LIS_logunit, *)'[INFO] AC72 Irrigation_Filename: not defined, no irrigation'
            AC72_struc(n)%Irrigation_Filename = '(None)'
        endif
    enddo
 
    ! AquaCrop model soil parameter table
    call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 soil parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC72_struc(n)%soil_tbl_name, rc=rc)
        call LIS_verify(rc, "AquaCrop.7.2 soil parameter table: not defined")
    enddo

 
    ! soil classification scheme
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'SOILTEXT_SCHEME', soil_scheme_name)
        call LIS_verify(ios, 'Error in nf90_get_att: SOILTEXT_SCHEME')
        if (trim(soil_scheme_name) .eq. "STATSGO") then
          AC72_struc(n)%soil_scheme_name = "STAS"
        else
          write(LIS_logunit, *) "Fatal error: currently, only STATSGO soil scheme is supported by AquaCrop!"
          call LIS_endrun()
        endif 
    enddo

    ! Read Nrsoil layers from LDT
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'SOIL_LAYERS', AC72_struc(n)%NrSoilLayers)
        call LIS_verify(ios,'Error in nf90_get_att in ac72_readcrd')
        write(LIS_logunit,*)"[INFO] Number of soil layers from LDT: ",AC72_struc(n)%NrSoilLayers
        allocate(AC72_struc(n)%Thickness(AC72_struc(n)%NrSoilLayers))
        allocate(AC72_struc(n)%init_smc(AC72_struc(n)%NrSoilLayers))
        
        ! Get thickness of soil layers from LDT
        do i=1, AC72_struc(n)%NrSoilLayers
            write (str, '(i0)') i
            ios = nf90_get_att(nids(n), NF90_GLOBAL, "THICKNESS_LAYER_"//trim(str), &
                                AC72_struc(n)%Thickness(i))
            call LIS_verify(ios,'Error in nf90_get_att in ac72_readcrd')
            write(LIS_logunit,*)"[INFO] AC72 THICKNESS_LAYER_"//trim(str), &
                                AC72_struc(n)%Thickness(i)
        enddo
    enddo

    ! Reference year for climatology (stress functions)
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'AC_CLIM_REF_YEAR', AC72_struc(n)%tempcli_refyr)
        call LIS_verify(ios,'Error in nf90_get_att in ac72_readcrd')
        write(LIS_logunit,*)"[INFO] Reference year for climatology from LDT: ",AC72_struc(n)%tempcli_refyr
    enddo
 
    do n=1, LIS_rc%nnest
      AC72_struc(n)%dt = AC72_struc(n)%ts
    enddo 

    ! reference height of forcings (T and q)
    call ESMF_ConfigFindLabel(LIS_config, &
         "AquaCrop.7.2 reference height of T and q:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC72_struc(n)%refz_tq, rc=rc)
        call LIS_verify(rc, &
             "AquaCrop.7.2 reference height of T and q: "//&
             "not defined")
    enddo

    ! reference height of forcings (u and v)
    call ESMF_ConfigFindLabel(LIS_config, &
         "AquaCrop.7.2 reference height of u and v:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, AC72_struc(n)%refz_uv, rc=rc)
        call LIS_verify(rc, &
             "AquaCrop.7.2 reference height of u and v: "//&
             "not defined")
    enddo

    do n=1,LIS_rc%nnest
      ios = nf90_close(nids(n))
      call LIS_verify(ios,'Error in nf90_close in '//trim(LIS_rc%paramfile(n))//' in AC72_readcrd')
    enddo 

    ! The following lines hard code the LDT NetCDF variable names. 
    do n=1, LIS_rc%nnest
        AC72_struc(n)%LDT_ncvar_soiltype = 'AC72_SOILTYPE'
        AC72_struc(n)%LDT_ncvar_tmincli_monthly = 'AC_Tmin_clim'
        AC72_struc(n)%LDT_ncvar_tmaxcli_monthly = 'AC_Tmax_clim'
    enddo

    ! set default restart format to netcdf
    do n=1,LIS_rc%nnest
        AC72_struc(n)%rformat = "netcdf"
    enddo
    ! restart run, read restart file
    if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC72_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.2 restart file: not defined")
        enddo

        Call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC72_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.2 restart file format: not defined")
        enddo
        ! cold start run, read initial state variables
    else
        ! volumetric liquid soil moisture
        call ESMF_ConfigFindLabel(LIS_config, "AquaCrop.7.2 initial liquid soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, AC72_struc(n)%init_smc, rc=rc)
            call LIS_verify(rc, "AquaCrop.7.2 initial liquid soil moistures: not defined")
        enddo
    endif

    deallocate(nids)

    write(LIS_logunit, *) "[INFO] Finish reading LIS configuration file for AquaCrop.7.2 model"
end subroutine AC72_readcrd
