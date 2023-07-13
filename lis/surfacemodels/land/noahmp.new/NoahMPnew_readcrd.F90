!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!BOP
!
! !ROUTINE: NoahMPnew_readcrd
! \label{NoahMPnew\_readcrd}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  10/25/18 : Shugong Wang, Zhuo Wang, initial implementation for LIS 7 and NoahMP401
!  May 2023: Cenlin He; update to work with refactored NoahMP (v5.0 and newer)
!
! !INTERFACE:
subroutine NoahMPnew_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod,      only : LIS_rc , LIS_config
    use LIS_timeMgrMod,   only : LIS_parseTimeString
    use LIS_logMod,       only : LIS_logunit, LIS_verify, LIS_endrun
    use NoahMPnew_lsmMod, only : NoahmpNew_struc
    use netcdf
!
! !DESCRIPTION:
!
!  This routine reads the options specific to NoahMPnew model from
!  the LIS configuration file.
!
!EOP
    implicit none

    integer      :: rc 
    integer      :: n, i
    character*10 :: time 
    character*6  :: str_i
    integer      :: ios
    integer, allocatable :: nids(:)
    character*32 :: soil_scheme_name, landuse_scheme_name

    allocate(nids(LIS_rc%nnest))
 
    write(LIS_logunit,*) &
         "[INFO] Start reading LIS configuration file for Noah-MP.New (v5.0 or newer)"
    
    ! open NetCDF parameter file for reading global attributes 
    do n=1,LIS_rc%nnest
      ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE,ncid=nids(n))
      call LIS_verify(ios,'Error in nf90_open in '//trim(LIS_rc%paramfile(n))//' in NoahMPnew_readcrd')
    enddo 

    ! main Noah-MP model timestep
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "Noah-MP.New model timestep: not defined")
        call LIS_parseTimeString(time, NoahmpNew_struc(n)%ts)
    enddo

    ! Noah-MP soil process timestep
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New soil timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "Noah-MP.New soil timestep: not defined")
        call LIS_parseTimeString(time, NoahmpNew_struc(n)%ts_soil)
    enddo
    
    ! restart timestep
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, &
             "Noah-MP.New restart output interval: not defined")
        call LIS_parseTimeString(time, NoahmpNew_struc(n)%rstInterval)
    enddo
    
    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number of soil layers
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New number of soil layers:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%nsoil, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New number of soil layers: not defined")
    enddo
 
    ! allocate memory for sldpth using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(NoahmpNew_struc(n)%sldpth(NoahmpNew_struc(n)%nsoil))
    enddo
    ! allocate memory for init_smc using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(NoahmpNew_struc(n)%init_smc(NoahmpNew_struc(n)%nsoil))
    enddo
    ! allocate memory for init_tslb using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(NoahmpNew_struc(n)%init_tslb(NoahmpNew_struc(n)%nsoil))
    enddo

    ! maximum number of snow layers (e.g., 3)
    do n=1, LIS_rc%nnest
       NoahmpNew_struc(n)%nsnow = 3
    enddo

    ! thickness of atmospheric layers
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New reference height of temperature and humidity:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%dz8w, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New reference height of temperature and "//&
             "humidity: not defined")
    enddo
 
    ! thickness of soil layers
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New thickness of soil layers:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, NoahmpNew_struc(n)%nsoil
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%sldpth(i), rc=rc)
            call LIS_verify(rc, &
                 'Noah-MP.New thickness of soil layers: not defined')
        enddo
    enddo
 
    ! Landuse classification scheme
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'LANDCOVER_SCHEME', landuse_scheme_name)
        call LIS_verify(ios, 'Error in nf90_get_att: LANDCOVER_SCHEME')
        if (trim(landuse_scheme_name) .eq. "USGS") then
          NoahmpNew_struc(n)%landuse_scheme_name = "USGS"
        elseif (trim(landuse_scheme_name) .eq. "IGBPNCEP") then
          NoahmpNew_struc(n)%landuse_scheme_name = &
                               "MODIFIED_IGBP_MODIS_NOAH"
        elseif (trim(landuse_scheme_name) .eq. "NALCMS_SM_IGBPNCEP" ) then
          NoahmpNew_struc(n)%landuse_scheme_name = &
                               "MODIFIED_IGBP_MODIS_NOAH"
        elseif (trim(landuse_scheme_name) .eq. "UMD") then
          NoahmpNew_struc(n)%landuse_scheme_name = "UMD"
        else
         write(LIS_logunit,*) &
                         "[ERR] Currently, only USGS, IGBPNCEP, and UMD"
         write(LIS_logunit,*) "[ERR] are supported by Noah-MP-New LSM"
         write(LIS_logunit,*) "[ERR] program stopping ..."
         call LIS_endrun()
        endif
    enddo
 
    ! NoahMP.New parameter table (merged SOILPARM.TBL,GENPARM.TBL,MPTABLE.TBL)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%noahmp_tbl_name, rc=rc)
        call LIS_verify(rc, &
        "Noah-MP.New parameter table: not defined")
    enddo

    write(LIS_logunit,*) &
          "[INFO] Setting Noah-MP.New physics options:"
 
    ! dynamic vegetation
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New dynamic vegetation option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%dveg_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New dynamic vegetation option: not defined")
        write(LIS_logunit,33) "dynamic vegetation:", &
                               NoahmpNew_struc(n)%dveg_opt
    enddo

    ! canopy stomatal resistance (1->Ball-Berry; 2->Jarvis)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New canopy stomatal resistance option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%crs_opt, rc=rc)
        call LIS_verify(rc, &
          "Noah-MP.New canopy stomatal resistance option: not defined")
        write(LIS_logunit,33) "canopy stomatal resistance:", &
                               NoahmpNew_struc(n)%crs_opt
    enddo

    ! soil moisture factor for stomatal resistance(1->Noah;2->CLM;3->SSiB)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New soil moisture factor for stomatal resistance:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%btr_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.New soil moisture factor for stomatal resistance:"//&
           " not defined")
        write(LIS_logunit,33) "soil moisture factor for stomatal "//&
                              "resistance:",NoahmpNew_struc(n)%btr_opt
    enddo

    ! surface runoff (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS; 5->MMF; 6->VIC; 7->XAJ; 8->DynVIC)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New surface runoff option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%runsfc_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New surface runoff option: not defined")
        write(LIS_logunit,33) "surface runoff:", &
                               NoahmpNew_struc(n)%runsfc_opt
    enddo

    ! subsurface runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS; 5->MMF; 6->VIC; 7->XAJ; 8->DynVIC)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New subsurface runoff and groundwater option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%runsub_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New subsurface runoff and groundwater option: not defined")
        write(LIS_logunit,33) "subsurface runoff and groundwater:", &
                               NoahmpNew_struc(n)%runsub_opt
    enddo

    ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New surface layer drag coefficient option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%sfc_opt, rc=rc)
        call LIS_verify(rc, &
            "Noah-MP.New surface layer drag coefficient option:"//&
           " not defined")
        write(LIS_logunit,33) "surface layer drag coefficient:", &
                              NoahmpNew_struc(n)%sfc_opt
    enddo

    ! supercooled liquid water (1->NY06; 2->Koren99)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New supercooled liquid water option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%frz_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.New supercooled liquid water option: not defined")
        write(LIS_logunit,33) "supercooled liquid water:", &
                              NoahmpNew_struc(n)%frz_opt
    enddo

    ! frozen soil permeability (1->NY06; 2->Koren99)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New frozen soil permeability option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%inf_opt, rc=rc)
        call LIS_verify(rc, &
            "Noah-MP.New frozen soil permeability option: not defined")
        write(LIS_logunit,33) "frozen soil permeability:", &
                              NoahmpNew_struc(n)%inf_opt
    enddo
 
    ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New canopy radiative transfer option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%rad_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New canopy radiative transfer option: not defined")
        write(LIS_logunit,33) "canopy radiative transfer:", &
                              NoahmpNew_struc(n)%rad_opt
    enddo
 
    ! snow surface albedo (1->BATS; 2->CLASS)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New snow surface albedo option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%alb_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New snow surface albedo option: not defined")
        write(LIS_logunit,33) "snow surface albedo:", &
                              NoahmpNew_struc(n)%alb_opt
    enddo
 
    ! rainfall & snowfall (1->Jordan91; 2->BATS; 3->Noah; 4->WRF; 5->Wet-bulb)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New rain-snow partition option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%snf_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New rain-snow partition option: not defined")
        write(LIS_logunit,33) "rain-snow partition:", &
                              NoahmpNew_struc(n)%snf_opt
    enddo
 
    ! lower boundary of soil temperature
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New lower boundary of soil temperature option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%tbot_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.New lower boundary of soil temperature option:"//&
           " not defined")
        write(LIS_logunit,33) "lower boundary of soil temperature:", &
                              NoahmpNew_struc(n)%tbot_opt
    enddo
 
    ! snow/soil temperature time scheme
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New snow&soil temperature time scheme option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%stc_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.New snow&soil temperature time scheme option:"//&
           " not defined")
        write(LIS_logunit,33) "snow&soil temperature time scheme:", &
                              NoahmpNew_struc(n)%stc_opt
    enddo
 
    ! glacier option (1->phase change; 2->simple)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New glacier ice option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%gla_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New glacier ice option: not defined")
        write(LIS_logunit,33) "glacier ice:",NoahmpNew_struc(n)%gla_opt
    enddo

    ! Custom snowpack depth for glacier model (in mm)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New snow depth glacier model option:", rc = rc)
    if(rc /= 0) then
        write(LIS_logunit,33) "[WARN] Max snow depth not defined."
        write(LIS_logunit,33) "[WARN] Setting to default value of 2000."
        do n=1, LIS_rc%nnest
            NoahmpNew_struc(n)%sndpth_gla_opt = 2000
            write(LIS_logunit,33) "snow depth for glacier model: ",NoahmpNew_struc(n)%sndpth_gla_opt
        enddo
    else
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%sndpth_gla_opt, rc=rc)
            write(LIS_logunit,33) "snow depth for glacier model: ",NoahmpNew_struc(n)%sndpth_gla_opt
        enddo
    endif

    ! surface resistance (1->Sakaguchi/Zeng;2->Seller;3->mod Sellers;4->1+snow)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New surface resistance option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%rsf_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New surface resistance option: not defined")
        write(LIS_logunit,33) "surface resistance:", &
                              NoahmpNew_struc(n)%rsf_opt
    enddo
 
    ! soil configuration option
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New soil configuration option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%soil_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New soil configuration option: not defined")
        write(LIS_logunit,33) "soil configuration:", &
                              NoahmpNew_struc(n)%soil_opt
    enddo
 
    ! soil pedotransfer function option
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New soil pedotransfer function option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%pedo_opt, rc=rc)
        call LIS_verify(rc, &
          "Noah-MP.New soil pedotransfer function option: not defined")
        write(LIS_logunit,33) "soil pedotransfer function:", &
                              NoahmpNew_struc(n)%pedo_opt
    enddo
 
    ! crop model option (0->none; 1->Liu2016)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New crop model option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%crop_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New crop model option: not defined")
        write(LIS_logunit,33) "crop model:", &
                              NoahmpNew_struc(n)%crop_opt
    enddo
 
    ! urban physics option
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New urban physics option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%urban_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New urban physics option: not defined")
        write(LIS_logunit,33) "urban physics:", &
                              NoahmpNew_struc(n)%urban_opt
    enddo

    ! snow thermal conductivity option (1->Yen1965; 2->Anderson1976; 3->Constant; 4->Verseghy1991; 5->Yen1981)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New snow thermal conductivity option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%tksno_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New snow thermal conductivity option: not defined")
        write(LIS_logunit,33) "snow thermal conductivity:",NoahmpNew_struc(n)%tksno_opt
    enddo

    ! irrigation option (0->none; 1->always on; 2->trigger by planting/harvest dates; 3->trigger by LAI)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New irrigation trigger option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%irr_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New irrigation trigger option: not defined")
        write(LIS_logunit,33) "irrigation trigger:",NoahmpNew_struc(n)%irr_opt
    enddo

    ! irrigation method option (0->fraction from input; 1->sprinkler; 2->micro/drip; 3->flood)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New irrigation method option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%irrm_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New irrigation method option: not defined")
        write(LIS_logunit,33) "irrigation method:",NoahmpNew_struc(n)%irrm_opt
    enddo

    ! tile drainage option (0->none; 1->simple drainage; 2->Hooghoudt's scheme)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.New tile drainage option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%tdrn_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.New tile drainage option: not defined")
        write(LIS_logunit,33) "tile drainage:",NoahmpNew_struc(n)%tdrn_opt
    enddo


    ! The following lines hard code the LDT NetCDF variable names. 
    ! Modified by Zhuo Wang on 11/08/2018
    ! Setting some values to PLANTING, HARVEST, SEASON_GDD, SOILCOMP, SOILCL1-->SOILCL4 in NoahMPnew_main.F90
    do n=1, LIS_rc%nnest
        ! NoahmpNew_struc(n)%LDT_ncvar_vegetype = ' ! Edit here if hard code name
        ! NoahmpNew_struc(n)%LDT_ncvar_soiltype = ' ! Edit here if hard code name
        NoahmpNew_struc(n)%LDT_ncvar_tbot    = 'TBOT'               !'NOAHMPnew_TBOT'
        NoahmpNew_struc(n)%LDT_ncvar_shdfac_monthly = 'GREENNESS'   !'NOAHMPnew_SHDFAC_MONTHLY'
        NoahmpNew_struc(n)%LDT_ncvar_planting = 'PLANTING'          !'NOAHMPnew_PLANTING'
        NoahmpNew_struc(n)%LDT_ncvar_harvest = 'HARVEST'            !'NOAHMPnew_HARVEST'
        NoahmpNew_struc(n)%LDT_ncvar_season_gdd = 'SEASON_GDD'      !'NOAHMPnew_SEASON_GDD'
        NoahmpNew_struc(n)%LDT_ncvar_soilcomp = 'SOILCOMP'          !'NOAHMPnew_SOILCOMP'
        NoahmpNew_struc(n)%LDT_ncvar_soilcL1 = 'SOILCL1'            !'NOAHMPnew_SOILCL1'
        NoahmpNew_struc(n)%LDT_ncvar_soilcL2 = 'SOILCL2'            !'NOAHMPnew_SOILCL2'
        NoahmpNew_struc(n)%LDT_ncvar_soilcL3 = 'SOILCL3'            !'NOAHMPnew_SOILCL3'
        NoahmpNew_struc(n)%LDT_ncvar_soilcL4 = 'SOILCL4'            !'NOAHMPnew_SOILCL4'
    enddo

!------------------------------------------------------------------------------------------
    ! set default restart format to netcdf
    do n=1,LIS_rc%nnest
        NoahmpNew_struc(n)%rformat = "netcdf"
    enddo
    ! restart run, read restart file
    if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.New restart file: not defined")
        enddo
        
        Call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.New restart file format: not defined")
        enddo

    ! coldstart run, read initial state variables
    else
        ! skin temperature
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New initial surface skin temperature:", rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%init_tskin, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.New initial surface skin temperature:"//&
                 " not defined")
        enddo

        ! snow water equivalent
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New initial snow water equivalent:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%init_sneqv, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.New initial snow water equivalent:"//&
                 " not defined")
        enddo

        ! physical snow depth
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New initial snow depth:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%init_snowh, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.New initial snow depth:"//&
                 " not defined")
        enddo

        ! total canopy water + ice
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New initial total canopy surface water:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%init_canwat, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.New initial total canopy surface water:"//&
                 " not defined")
        enddo

        ! soil temperature
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New initial soil temperatures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, NoahmpNew_struc(n)%nsoil
               call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%init_tslb(i), rc=rc)
            end do
            call LIS_verify(rc, &
                 "Noah-MP.New initial soil temperatures:"//&
                 " not defined")
        enddo

        ! volumetric soil moisture
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New initial total soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, NoahmpNew_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%init_smc(i), rc=rc)
            end do
            call LIS_verify(rc, &
                 "Noah-MP.New initial total soil moistures:"//&
                 " not defined")
        enddo

        ! water table depth
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New initial water table depth:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%init_zwt, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.New initial water table depth:"//&
                 " not defined")
        enddo

        ! water in the "aquifer"
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New initial water in the aquifer:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%init_wa, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.New initial water in the aquifer:"//&
                 " not defined")
        enddo

        ! water in aquifer and saturated soil
        call ESMF_ConfigFindLabel(LIS_config, &
            "Noah-MP.New initial water in aquifer and saturated soil:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%init_wt, rc=rc)
            call LIS_verify(rc, &
            "Noah-MP.New initial water in aquifer and saturated soil:"//&
                 " not defined")
        enddo

        ! leaf area index
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.New initial leaf area index:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NoahmpNew_struc(n)%init_lai, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.New initial leaf area index:"//&
                 " not defined")
        enddo

    endif
     
    deallocate(nids)

 33   format(a47,i4)

    write(LIS_logunit, *) &
        "[INFO] Finish reading LIS configuration file for Noah-MP.New"
     
end subroutine NoahMPnew_readcrd
