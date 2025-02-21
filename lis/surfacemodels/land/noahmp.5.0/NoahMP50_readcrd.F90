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
! !ROUTINE: NoahMP50_readcrd
! \label{NoahMP50\_readcrd}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!  10/25/18 : Shugong Wang, Zhuo Wang, initial implementation for LIS 7 and NoahMP401
!  May 2023: Cenlin He; update to work with refactored NoahMP (v5.0 and newer)
!
! !INTERFACE:

subroutine NoahMP50_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod,      only : LIS_rc , LIS_config
    use LIS_timeMgrMod,   only : LIS_parseTimeString
    use LIS_logMod,       only : LIS_logunit, LIS_verify, LIS_endrun
    use NoahMP50_lsmMod, only : Noahmp50_struc
    use netcdf
!
! !DESCRIPTION:
!
!  This routine reads the options specific to NoahMP50 model from
!  the LIS configuration file.
!
!EOP
    implicit none

    integer      :: rc 
    integer      :: n, i
    character*10 :: time 
    integer      :: ios
    integer, allocatable :: nids(:)
    character*32 :: landuse_scheme_name

    allocate(nids(LIS_rc%nnest))
 
    write(LIS_logunit,*) &
         "[INFO] Start reading LIS configuration file for Noah-MP.5.0 (v5.0 or newer)"
    
    ! open NetCDF parameter file for reading global attributes 
    do n=1,LIS_rc%nnest
      ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE,ncid=nids(n))
      call LIS_verify(ios,'Error in nf90_open in '//trim(LIS_rc%paramfile(n))//' in NoahMP50_readcrd')
    enddo 

    ! main Noah-MP model timestep
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "Noah-MP.5.0 model timestep: not defined")
        call LIS_parseTimeString(time, Noahmp50_struc(n)%ts)
    enddo

    ! Noah-MP soil process timestep
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 soil timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "Noah-MP.5.0 soil timestep: not defined")
        call LIS_parseTimeString(time, Noahmp50_struc(n)%ts_soil)
    enddo
    
    ! restart timestep
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 restart output interval: not defined")
        call LIS_parseTimeString(time, Noahmp50_struc(n)%rstInterval)
    enddo

    ! model domain size dx (meter)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 domain resolution dx:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%dx, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 domain resolution dx: not defined")
    enddo

    ! model domain size dy (meter)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 domain resolution dy:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%dy, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 domain resolution dy: not defined")
    enddo

    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number of soil layers
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 number of soil layers:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%nsoil, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 number of soil layers: not defined")
    enddo
 
    ! allocate memory for sldpth using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(Noahmp50_struc(n)%sldpth(Noahmp50_struc(n)%nsoil))
    enddo
    ! allocate memory for init_smc using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(Noahmp50_struc(n)%init_smc(Noahmp50_struc(n)%nsoil))
    enddo
    ! allocate memory for init_tslb using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(Noahmp50_struc(n)%init_tslb(Noahmp50_struc(n)%nsoil))
    enddo

    ! maximum number of snow layers (e.g., 3)
    do n=1, LIS_rc%nnest
       Noahmp50_struc(n)%nsnow = 3
    enddo

    ! thickness of atmospheric layers
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 reference height of temperature and humidity:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%dz8w, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 reference height of temperature and "//&
             "humidity: not defined")
    enddo
 
    ! thickness of soil layers
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 thickness of soil layers:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, Noahmp50_struc(n)%nsoil
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%sldpth(i), rc=rc)
            call LIS_verify(rc, &
                 'Noah-MP.5.0 thickness of soil layers: not defined')
        enddo
    enddo
 
    ! Landuse classification scheme
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'LANDCOVER_SCHEME', landuse_scheme_name)
        call LIS_verify(ios, 'Error in nf90_get_att: LANDCOVER_SCHEME')
        if (trim(landuse_scheme_name) .eq. "USGS") then
          Noahmp50_struc(n)%landuse_scheme_name = "USGS"
        elseif (trim(landuse_scheme_name) .eq. "IGBPNCEP") then
          Noahmp50_struc(n)%landuse_scheme_name = &
                               "MODIFIED_IGBP_MODIS_NOAH"
        elseif (trim(landuse_scheme_name) .eq. "NALCMS_SM_IGBPNCEP" ) then
          Noahmp50_struc(n)%landuse_scheme_name = &
                               "MODIFIED_IGBP_MODIS_NOAH"
        elseif (trim(landuse_scheme_name) .eq. "UMD") then
          Noahmp50_struc(n)%landuse_scheme_name = "UMD"
        else
         write(LIS_logunit,*) &
                         "[ERR] Currently, only USGS, IGBPNCEP, and UMD"
         write(LIS_logunit,*) "[ERR] are supported by Noah-MP.5.0 LSM"
         write(LIS_logunit,*) "[ERR] program stopping ..."
         call LIS_endrun()
        endif
    enddo
 
    ! Noah-MP.5.0 parameter table (merged SOILPARM.TBL,GENPARM.TBL,MPTABLE.TBL)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%noahmp_tbl_name, rc=rc)
        call LIS_verify(rc, &
        "Noah-MP.5.0 parameter table: not defined")
    enddo

    write(LIS_logunit,*) &
          "[INFO] Setting Noah-MP.5.0 physics options:"
 
    ! dynamic vegetation
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 dynamic vegetation option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%dveg_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 dynamic vegetation option: not defined")
        write(LIS_logunit,33) "dynamic vegetation:", &
                               Noahmp50_struc(n)%dveg_opt
    enddo

    ! canopy stomatal resistance (1->Ball-Berry; 2->Jarvis)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 canopy stomatal resistance option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%crs_opt, rc=rc)
        call LIS_verify(rc, &
          "Noah-MP.5.0 canopy stomatal resistance option: not defined")
        write(LIS_logunit,33) "canopy stomatal resistance:", &
                               Noahmp50_struc(n)%crs_opt
    enddo

    ! soil moisture factor for stomatal resistance(1->Noah;2->CLM;3->SSiB)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 soil moisture factor for stomatal resistance:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%btr_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.5.0 soil moisture factor for stomatal resistance:"//&
           " not defined")
        write(LIS_logunit,33) "soil moisture factor for stomatal "//&
                              "resistance:",Noahmp50_struc(n)%btr_opt
    enddo

    ! surface runoff (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS; 5->MMF; 6->VIC; 7->XAJ; 8->DynVIC)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 surface runoff option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%runsfc_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 surface runoff option: not defined")
        write(LIS_logunit,33) "surface runoff:", &
                               Noahmp50_struc(n)%runsfc_opt
    enddo

    ! subsurface runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS; 5->MMF; 6->VIC; 7->XAJ; 8->DynVIC)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 subsurface runoff and groundwater option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%runsub_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 subsurface runoff and groundwater option: not defined")
        write(LIS_logunit,33) "subsurface runoff and groundwater:", &
                               Noahmp50_struc(n)%runsub_opt
    enddo

    ! infiltration options for dynamic VIC (1->Philip; 2-> Green-Ampt;3->Smith-Parlange)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 dynamic VIC infiltration option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%infdv_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 dynamic VIC infiltration option: not defined")
        write(LIS_logunit,33) "dynamic VIC infiltration:", &
                               Noahmp50_struc(n)%infdv_opt
    enddo

    ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 surface layer drag coefficient option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%sfc_opt, rc=rc)
        call LIS_verify(rc, &
            "Noah-MP.5.0 surface layer drag coefficient option:"//&
           " not defined")
        write(LIS_logunit,33) "surface layer drag coefficient:", &
                              Noahmp50_struc(n)%sfc_opt
    enddo

    ! supercooled liquid water (1->NY06; 2->Koren99)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 supercooled liquid water option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%frz_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.5.0 supercooled liquid water option: not defined")
        write(LIS_logunit,33) "supercooled liquid water:", &
                              Noahmp50_struc(n)%frz_opt
    enddo

    ! frozen soil permeability (1->NY06; 2->Koren99)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 frozen soil permeability option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%inf_opt, rc=rc)
        call LIS_verify(rc, &
            "Noah-MP.5.0 frozen soil permeability option: not defined")
        write(LIS_logunit,33) "frozen soil permeability:", &
                              Noahmp50_struc(n)%inf_opt
    enddo
 
    ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 canopy radiative transfer option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%rad_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 canopy radiative transfer option: not defined")
        write(LIS_logunit,33) "canopy radiative transfer:", &
                              Noahmp50_struc(n)%rad_opt
    enddo
 
    ! snow surface albedo (1->BATS; 2->CLASS)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 snow surface albedo option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%alb_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 snow surface albedo option: not defined")
        write(LIS_logunit,33) "snow surface albedo:", &
                              Noahmp50_struc(n)%alb_opt
    enddo
 
    ! rainfall & snowfall (1->Jordan91; 2->BATS; 3->Noah; 4->WRF; 5->Wet-bulb)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 rain-snow partition option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%snf_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 rain-snow partition option: not defined")
        write(LIS_logunit,33) "rain-snow partition:", &
                              Noahmp50_struc(n)%snf_opt
    enddo
 
    ! lower boundary of soil temperature
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 lower boundary of soil temperature option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%tbot_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.5.0 lower boundary of soil temperature option:"//&
           " not defined")
        write(LIS_logunit,33) "lower boundary of soil temperature:", &
                              Noahmp50_struc(n)%tbot_opt
    enddo
 
    ! snow/soil temperature time scheme
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 snow&soil temperature time scheme option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%stc_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.5.0 snow&soil temperature time scheme option:"//&
           " not defined")
        write(LIS_logunit,33) "snow&soil temperature time scheme:", &
                              Noahmp50_struc(n)%stc_opt
    enddo
 
    ! glacier option (1->phase change; 2->simple)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 glacier ice option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%gla_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 glacier ice option: not defined")
        write(LIS_logunit,33) "glacier ice:",Noahmp50_struc(n)%gla_opt
    enddo

    ! Custom snowpack depth for glacier model (in mm)
    !=== Now this parameter is set in the NoahmpTable.TBL ("SWEMAXGLA", default=5000mm)
    !call ESMF_ConfigFindLabel(LIS_config, &
    !     "Noah-MP.5.0 snow depth glacier model option:", rc = rc)
    !if(rc /= 0) then
    !    write(LIS_logunit,33) "[WARN] Max snow depth not defined."
    !    write(LIS_logunit,33) "[WARN] Setting to default value of 5000."
    !    do n=1, LIS_rc%nnest
    !        Noahmp50_struc(n)%sndpth_gla_opt = 5000
    !        write(LIS_logunit,33) "snow depth for glacier model: ",Noahmp50_struc(n)%sndpth_gla_opt
    !    enddo
    !else
    !    do n=1, LIS_rc%nnest
    !        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%sndpth_gla_opt, rc=rc)
    !        write(LIS_logunit,33) "snow depth for glacier model: ",Noahmp50_struc(n)%sndpth_gla_opt
    !    enddo
    !endif

    ! surface resistance (1->Sakaguchi/Zeng;2->Seller;3->mod Sellers;4->1+snow)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 surface resistance option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%rsf_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 surface resistance option: not defined")
        write(LIS_logunit,33) "surface resistance:", &
                              Noahmp50_struc(n)%rsf_opt
    enddo
 
    ! soil configuration option
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 soil configuration option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%soil_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 soil configuration option: not defined")
        write(LIS_logunit,33) "soil configuration:", &
                              Noahmp50_struc(n)%soil_opt
    enddo
 
    ! soil pedotransfer function option
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 soil pedotransfer function option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%pedo_opt, rc=rc)
        call LIS_verify(rc, &
          "Noah-MP.5.0 soil pedotransfer function option: not defined")
        write(LIS_logunit,33) "soil pedotransfer function:", &
                              Noahmp50_struc(n)%pedo_opt
    enddo
 
    ! crop model option (0->none; 1->Liu2016)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 crop model option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%crop_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 crop model option: not defined")
        write(LIS_logunit,33) "crop model:", &
                              Noahmp50_struc(n)%crop_opt
    enddo
 
    ! urban physics option
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 urban physics option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%urban_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 urban physics option: not defined")
        write(LIS_logunit,33) "urban physics:", &
                              Noahmp50_struc(n)%urban_opt
    enddo

    ! snow thermal conductivity option (1->Yen1965; 2->Anderson1976; 3->Constant; 4->Verseghy1991; 5->Yen1981)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 snow thermal conductivity option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%tksno_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 snow thermal conductivity option: not defined")
        write(LIS_logunit,33) "snow thermal conductivity:",Noahmp50_struc(n)%tksno_opt
    enddo

    ! irrigation option (0->none; 1->always on; 2->trigger by planting/harvest dates; 3->trigger by LAI)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 irrigation trigger option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%irr_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 irrigation trigger option: not defined")
        write(LIS_logunit,33) "irrigation trigger:",Noahmp50_struc(n)%irr_opt
    enddo

    ! irrigation method option (0->fraction from input; 1->sprinkler; 2->micro/drip; 3->flood)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 irrigation method option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%irrm_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 irrigation method option: not defined")
        write(LIS_logunit,33) "irrigation method:",Noahmp50_struc(n)%irrm_opt
    enddo

    ! tile drainage option (0->none; 1->simple drainage; 2->Hooghoudt's scheme)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.5.0 tile drainage option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%tdrn_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.5.0 tile drainage option: not defined")
        write(LIS_logunit,33) "tile drainage:",Noahmp50_struc(n)%tdrn_opt
    enddo


    ! The following lines hard code the LDT NetCDF variable names. 
    ! Setting some values to PLANTING, HARVEST, SEASON_GDD, SOILCOMP, SOILCL1-->SOILCL4 in NoahMP50_main.F90
    do n=1, LIS_rc%nnest
        ! Noahmp50_struc(n)%LDT_ncvar_vegetype = ' ! Edit here if hard code name
        ! Noahmp50_struc(n)%LDT_ncvar_soiltype = ' ! Edit here if hard code name
        Noahmp50_struc(n)%LDT_ncvar_tbot           = 'TBOT'             !'NoahMP50_TBOT'
        Noahmp50_struc(n)%LDT_ncvar_shdfac_monthly = 'GREENNESS'        !'NoahMP50_SHDFAC_MONTHLY'
        Noahmp50_struc(n)%LDT_ncvar_planting       = 'PLANTING'         !'NoahMP50_PLANTING'
        Noahmp50_struc(n)%LDT_ncvar_harvest        = 'HARVEST'          !'NoahMP50_HARVEST'
        Noahmp50_struc(n)%LDT_ncvar_season_gdd     = 'SEASON_GDD'       !'NoahMP50_SEASON_GDD'
        Noahmp50_struc(n)%LDT_ncvar_soilcomp       = 'SOILCOMP'         !'NoahMP50_SOILCOMP'
        Noahmp50_struc(n)%LDT_ncvar_soilcL1        = 'SOILCL1'          !'NoahMP50_SOILCL1'
        Noahmp50_struc(n)%LDT_ncvar_soilcL2        = 'SOILCL2'          !'NoahMP50_SOILCL2'
        Noahmp50_struc(n)%LDT_ncvar_soilcL3        = 'SOILCL3'          !'NoahMP50_SOILCL3'
        Noahmp50_struc(n)%LDT_ncvar_soilcL4        = 'SOILCL4'          !'NoahMP50_SOILCL4'
        Noahmp50_struc(n)%LDT_ncvar_irfract        = 'IRFRACT'
        Noahmp50_struc(n)%LDT_ncvar_sifract        = 'SIFRACT'
        Noahmp50_struc(n)%LDT_ncvar_mifract        = 'MIFRACT'
        Noahmp50_struc(n)%LDT_ncvar_fifract        = 'FIFRACT'
        Noahmp50_struc(n)%LDT_ncvar_tdfract        = 'TD_FRACTION'
        Noahmp50_struc(n)%LDT_ncvar_fdepth         = 'FDEPTH'
        Noahmp50_struc(n)%LDT_ncvar_eqzwt          = 'EQZWT'
        Noahmp50_struc(n)%LDT_ncvar_rechclim       = 'RECHCLIM'
        Noahmp50_struc(n)%LDT_ncvar_riverbed       = 'RIVERBED'
    enddo

!------------------------------------------------------------------------------------------
    ! set default restart format to netcdf
    do n=1,LIS_rc%nnest
        Noahmp50_struc(n)%rformat = "netcdf"
    enddo
    ! restart run, read restart file
    if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.5.0 restart file: not defined")
        enddo
        
        Call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.5.0 restart file format: not defined")
        enddo

    ! coldstart run, read initial state variables
    else
        ! skin temperature
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 initial surface skin temperature:", rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%init_tskin, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.5.0 initial surface skin temperature:"//&
                 " not defined")
        enddo

        ! snow water equivalent
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 initial snow water equivalent:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%init_sneqv, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.5.0 initial snow water equivalent:"//&
                 " not defined")
        enddo

        ! physical snow depth
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 initial snow depth:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%init_snowh, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.5.0 initial snow depth:"//&
                 " not defined")
        enddo

        ! total canopy water + ice
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 initial total canopy surface water:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%init_canwat, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.5.0 initial total canopy surface water:"//&
                 " not defined")
        enddo

        ! soil temperature
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 initial soil temperatures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, Noahmp50_struc(n)%nsoil
               call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%init_tslb(i), rc=rc)
            end do
            call LIS_verify(rc, &
                 "Noah-MP.5.0 initial soil temperatures:"//&
                 " not defined")
        enddo

        ! volumetric soil moisture
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 initial total soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, Noahmp50_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%init_smc(i), rc=rc)
            end do
            call LIS_verify(rc, &
                 "Noah-MP.5.0 initial total soil moistures:"//&
                 " not defined")
        enddo

        ! water table depth
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 initial water table depth:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%init_zwt, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.5.0 initial water table depth:"//&
                 " not defined")
        enddo

        ! water in the "aquifer"
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 initial water in the aquifer:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%init_wa, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.5.0 initial water in the aquifer:"//&
                 " not defined")
        enddo

        ! water in aquifer and saturated soil
        call ESMF_ConfigFindLabel(LIS_config, &
            "Noah-MP.5.0 initial water in aquifer and saturated soil:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%init_wt, rc=rc)
            call LIS_verify(rc, &
            "Noah-MP.5.0 initial water in aquifer and saturated soil:"//&
                 " not defined")
        enddo

        ! leaf area index
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.5.0 initial leaf area index:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, Noahmp50_struc(n)%init_lai, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.5.0 initial leaf area index:"//&
                 " not defined")
        enddo

    endif
     
    deallocate(nids)

 33   format(a47,i4)

    write(LIS_logunit, *) &
        "[INFO] Finish reading LIS configuration file for Noah-MP.5.0"
     
end subroutine NoahMP50_readcrd
