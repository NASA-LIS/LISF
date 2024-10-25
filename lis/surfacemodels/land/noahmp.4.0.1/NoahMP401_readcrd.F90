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
! !ROUTINE: NoahMP401_readcrd
! \label{NoahMP401\_readcrd}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/25/18 : Shugong Wang, Zhuo Wang, initial implementation for LIS 7 and NoahMP401
!
! !INTERFACE:
subroutine NoahMP401_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit, LIS_verify, LIS_endrun
    use NoahMP401_lsmMod, only       : NOAHMP401_struc
    use netcdf
!
! !DESCRIPTION:
!
!  This routine reads the options specific to NoahMP401 model from
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
         "[INFO] Start reading LIS configuration file for Noah-MP.4.0.1"
    
    ! open NetCDF parameter file for reading global attributes 
    do n=1,LIS_rc%nnest
      ios = nf90_open(path=trim(LIS_rc%paramfile(n)), mode=NF90_NOWRITE,ncid=nids(n))
      call LIS_verify(ios,'Error in nf90_open in '//trim(LIS_rc%paramfile(n))//' in NoahMP401_readcrd')
    enddo 

    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "Noah-MP.4.0.1 model timestep: not defined")
        call LIS_parseTimeString(time, NOAHMP401_struc(n)%ts)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 restart output interval: not defined")
        call LIS_parseTimeString(time, NOAHMP401_struc(n)%rstInterval)
    enddo
    
    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number of soil layers
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 number of soil layers:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%nsoil, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 number of soil layers: not defined")
    enddo
 
    ! allocate memory for sldpth using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(NOAHMP401_struc(n)%sldpth(NOAHMP401_struc(n)%nsoil))
    enddo
    ! allocate memory for init_smc using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(NOAHMP401_struc(n)%init_smc(NOAHMP401_struc(n)%nsoil))
    enddo
    ! allocate memory for init_tslb using nsoil as dimension
    do n=1, LIS_rc%nnest
        allocate(NOAHMP401_struc(n)%init_tslb(NOAHMP401_struc(n)%nsoil))
    enddo

    ! maximum number of snow layers (e.g., 3)
    do n=1, LIS_rc%nnest
       NOAHMP401_struc(n)%nsnow = 3
    enddo

    ! allocate memory for init_gecros_state with a dimension of 60
    do n=1, LIS_rc%nnest
        allocate(NOAHMP401_struc(n)%init_gecros_state(60))
    enddo
 
    ! thickness of atmospheric layers
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 reference height of temperature and humidity:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%dz8w, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 reference height of temperature and "//&
             "humidity: not defined")
    enddo
 
    ! thickness of soil layers
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 thickness of soil layers:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, NOAHMP401_struc(n)%nsoil
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%sldpth(i), rc=rc)
            call LIS_verify(rc, &
                 'Noah-MP.4.0.1 thickness of soil layers: not defined')
        enddo
    enddo
 
    ! Landuse classification scheme
    do n=1, LIS_rc%nnest
        ios = nf90_get_att(nids(n), NF90_GLOBAL, 'LANDCOVER_SCHEME', landuse_scheme_name)
        call LIS_verify(ios, 'Error in nf90_get_att: LANDCOVER_SCHEME')
        if (trim(landuse_scheme_name) .eq. "USGS") then
          NOAHMP401_struc(n)%landuse_scheme_name = "USGS"
        elseif (trim(landuse_scheme_name) .eq. "IGBPNCEP") then
          NOAHMP401_struc(n)%landuse_scheme_name = &
                               "MODIFIED_IGBP_MODIS_NOAH"
        elseif (trim(landuse_scheme_name) .eq. "NALCMS_SM_IGBPNCEP" ) then
          NOAHMP401_struc(n)%landuse_scheme_name = &
                               "MODIFIED_IGBP_MODIS_NOAH"
        elseif (trim(landuse_scheme_name) .eq. "UMD") then
          NOAHMP401_struc(n)%landuse_scheme_name = "UMD"
        else
         write(LIS_logunit,*) &
                         "[ERR] Currently, only USGS, IGBPNCEP, and UMD"
         write(LIS_logunit,*) "[ERR] are supported by Noah-MP-4.0.1 LSM"
         write(LIS_logunit,*) "[ERR] program stopping ..."
         call LIS_endrun()
        endif
    enddo
 
    ! Noah model soil parameter table
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 soil parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%soil_tbl_name, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 soil parameter table: not defined")
    enddo
 
    ! Noah model general parameter table
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 general parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%gen_tbl_name, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 general parameter table: not defined")
    enddo
 
    ! NoahMP parameter table
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 MP parameter table:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%noahmp_tbl_name, rc=rc)
        call LIS_verify(rc, &
        "Noah-MP.4.0.1 MP parameter table: not defined")
    enddo

    write(LIS_logunit,*) &
          "[INFO] Setting Noah-MP.4.0.1 physics options:"
 
    ! dynamic vegetation
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 dynamic vegetation option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%dveg_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 dynamic vegetation option: not defined")
        write(LIS_logunit,33) "dynamic vegetation:", &
                               NOAHMP401_struc(n)%dveg_opt
    enddo

    ! canopy stomatal resistance (1->Ball-Berry; 2->Jarvis)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 canopy stomatal resistance option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%crs_opt, rc=rc)
        call LIS_verify(rc, &
          "Noah-MP.4.0.1 canopy stomatal resistance option: not defined")
        write(LIS_logunit,33) "canopy stomatal resistance:", &
                               NOAHMP401_struc(n)%crs_opt
    enddo

    ! soil moisture factor for stomatal resistance(1->Noah;2->CLM;3->SSiB)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 soil moisture factor for stomatal resistance:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%btr_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.4.0.1 soil moisture factor for stomatal resistance:"//&
           " not defined")
        write(LIS_logunit,33) "soil moisture factor for stomatal "//&
                              "resistance:",NOAHMP401_struc(n)%btr_opt
    enddo

    ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 runoff and groundwater option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%run_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 runoff and groundwater option: not defined")
        write(LIS_logunit,33) "runoff and groundwater:", &
                               NOAHMP401_struc(n)%run_opt
    enddo

    ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 surface layer drag coefficient option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%sfc_opt, rc=rc)
        call LIS_verify(rc, &
            "Noah-MP.4.0.1 surface layer drag coefficient option:"//&
           " not defined")
        write(LIS_logunit,33) "surface layer drag coefficient:", &
                              NOAHMP401_struc(n)%sfc_opt
    enddo

    ! supercooled liquid water (1->NY06; 2->Koren99)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 supercooled liquid water option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%frz_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.4.0.1 supercooled liquid water option: not defined")
        write(LIS_logunit,33) "supercooled liquid water:", &
                              NOAHMP401_struc(n)%frz_opt
    enddo

    ! frozen soil permeability (1->NY06; 2->Koren99)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 frozen soil permeability option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%inf_opt, rc=rc)
        call LIS_verify(rc, &
            "Noah-MP.4.0.1 frozen soil permeability option: not defined")
        write(LIS_logunit,33) "frozen soil permeability:", &
                              NOAHMP401_struc(n)%inf_opt
    enddo
 
    ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 radiation transfer option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%rad_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 radiation transfer option: not defined")
        write(LIS_logunit,33) "radiation transfer:", &
                              NOAHMP401_struc(n)%rad_opt
    enddo
 
    ! snow surface albedo (1->BATS; 2->CLASS)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 snow surface albedo option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%alb_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 snow surface albedo option: not defined")
        write(LIS_logunit,33) "snow surface albedo:", &
                              NOAHMP401_struc(n)%alb_opt
    enddo
 
    ! rainfall & snowfall (1->Jordan91; 2->BATS; 3->Noah)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 rainfall & snowfall option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%snf_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 rainfall & snowfall option: not defined")
        write(LIS_logunit,33) "rainfall & snowfall:", &
                              NOAHMP401_struc(n)%snf_opt
    enddo
 
    ! lower boundary of soil temperature
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 lower boundary of soil temperature option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%tbot_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.4.0.1 lower boundary of soil temperature option:"//&
           " not defined")
        write(LIS_logunit,33) "lower boundary of soil temperature:", &
                              NOAHMP401_struc(n)%tbot_opt
    enddo
 
    ! snow/soil temperature time scheme
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 snow&soil temperature time scheme option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%stc_opt, rc=rc)
        call LIS_verify(rc, &
           "Noah-MP.4.0.1 snow&soil temperature time scheme option:"//&
           " not defined")
        write(LIS_logunit,33) "snow&soil temperature time scheme:", &
                              NOAHMP401_struc(n)%stc_opt
    enddo
 
    ! glacier option (1->phase change; 2->simple)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 glacier option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%gla_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 glacier option: not defined")
        write(LIS_logunit,33) "glacier:",NOAHMP401_struc(n)%gla_opt
    enddo

    ! Custom snowpack depth for glacier model (in mm)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 snow depth glacier model option:", rc = rc)
    if(rc /= 0) then
        write(LIS_logunit,33) "[WARN] Max snow depth not defined."
        write(LIS_logunit,33) "[WARN] Setting to default value of 2000."
        do n=1, LIS_rc%nnest
            NOAHMP401_struc(n)%sndpth_gla_opt = 2000
            write(LIS_logunit,33) "snow depth for glacier model: ",NOAHMP401_struc(n)%sndpth_gla_opt
        enddo
    else
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%sndpth_gla_opt, rc=rc)
            write(LIS_logunit,33) "snow depth for glacier model: ",NOAHMP401_struc(n)%sndpth_gla_opt
        enddo
    endif

    ! surface resistance (1->Sakaguchi/Zeng;2->Seller;3->mod Sellers;4->1+snow)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 surface resistance option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%rsf_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 surface resistance option: not defined")
        write(LIS_logunit,33) "surface resistance:", &
                              NOAHMP401_struc(n)%rsf_opt
    enddo
 
    ! soil configuration option
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 soil configuration option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%soil_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 soil configuration option: not defined")
        write(LIS_logunit,33) "soil configuration:", &
                              NOAHMP401_struc(n)%soil_opt
    enddo
 
    ! soil pedotransfer function option
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 soil pedotransfer function option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%pedo_opt, rc=rc)
        call LIS_verify(rc, &
          "Noah-MP.4.0.1 soil pedotransfer function option: not defined")
        write(LIS_logunit,33) "soil pedotransfer function:", &
                              NOAHMP401_struc(n)%pedo_opt
    enddo
 
    ! crop model option (0->none; 1->Liu et al.; 2->Gecros)
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 crop model option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%crop_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 crop model option: not defined")
        write(LIS_logunit,33) "crop model:", &
                              NOAHMP401_struc(n)%crop_opt
    enddo
 
    ! urban physics option
    call ESMF_ConfigFindLabel(LIS_config, &
         "Noah-MP.4.0.1 urban physics option:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%urban_opt, rc=rc)
        call LIS_verify(rc, &
             "Noah-MP.4.0.1 urban physics option: not defined")
        write(LIS_logunit,33) "urban physics:", &
                              NOAHMP401_struc(n)%urban_opt
    enddo

    ! The following lines hard code the LDT NetCDF variable names. 
    ! Modified by Zhuo Wang on 11/08/2018
    ! Setting some values to PLANTING, HARVEST, SEASON_GDD, SOILCOMP, SOILCL1-->SOILCL4 in NoahMP401_main.F90
    do n=1, LIS_rc%nnest
        ! NOAHMP401_struc(n)%LDT_ncvar_vegetype = ' ! Edit here if hard code name
        ! NOAHMP401_struc(n)%LDT_ncvar_soiltype = ' ! Edit here if hard code name
        NOAHMP401_struc(n)%LDT_ncvar_tbot    = 'TBOT'               !'NOAHMP401_TBOT'
        NOAHMP401_struc(n)%LDT_ncvar_shdfac_monthly = 'GREENNESS'   !'NOAHMP401_SHDFAC_MONTHLY'
        NOAHMP401_struc(n)%LDT_ncvar_planting = 'PLANTING'          !'NOAHMP401_PLANTING'
        NOAHMP401_struc(n)%LDT_ncvar_harvest = 'HARVEST'            !'NOAHMP401_HARVEST'
        NOAHMP401_struc(n)%LDT_ncvar_season_gdd = 'SEASON_GDD'      !'NOAHMP401_SEASON_GDD'
        NOAHMP401_struc(n)%LDT_ncvar_soilcomp = 'SOILCOMP'          !'NOAHMP401_SOILCOMP'
        NOAHMP401_struc(n)%LDT_ncvar_soilcL1 = 'SOILCL1'            !'NOAHMP401_SOILCL1'
        NOAHMP401_struc(n)%LDT_ncvar_soilcL2 = 'SOILCL2'            !'NOAHMP401_SOILCL2'
        NOAHMP401_struc(n)%LDT_ncvar_soilcL3 = 'SOILCL3'            !'NOAHMP401_SOILCL3'
        NOAHMP401_struc(n)%LDT_ncvar_soilcL4 = 'SOILCL4'            !'NOAHMP401_SOILCL4'
    enddo

!------------------------------------------------------------------------------------------
    ! set default restart format to netcdf
    do n=1,LIS_rc%nnest
        NOAHMP401_struc(n)%rformat = "netcdf"
    enddo
    ! restart run, read restart file
    if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 restart file: not defined")
        enddo
        
        Call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 restart file format: not defined")
        enddo

    ! coldstart run, read initial state variables
    else
        ! skin temperature
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 initial surface skin temperature:", rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_tskin, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 initial surface skin temperature:"//&
                 " not defined")
        enddo

        ! snow water equivalent
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 initial snow water equivalent:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_sneqv, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 initial snow water equivalent:"//&
                 " not defined")
        enddo

        ! physical snow depth
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 initial snow depth:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_snowh, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 initial snow depth:"//&
                 " not defined")
        enddo

        ! total canopy water + ice
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 initial total canopy surface water:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_canwat, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 initial total canopy surface water:"//&
                 " not defined")
        enddo

        ! soil temperature
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 initial soil temperatures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, NOAHMP401_struc(n)%nsoil
               call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_tslb(i), rc=rc)
            end do
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 initial soil temperatures:"//&
                 " not defined")
        enddo

        ! volumetric soil moisture
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 initial total soil moistures:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, NOAHMP401_struc(n)%nsoil
                call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_smc(i), rc=rc)
            end do
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 initial total soil moistures:"//&
                 " not defined")
        enddo

        ! water table depth
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 initial water table depth:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_zwt, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 initial water table depth:"//&
                 " not defined")
        enddo

        ! water in the "aquifer"
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 initial water in the aquifer:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_wa, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 initial water in the aquifer:"//&
                 " not defined")
        enddo

        ! water in aquifer and saturated soil
        call ESMF_ConfigFindLabel(LIS_config, &
            "Noah-MP.4.0.1 initial water in aquifer and saturated soil:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_wt, rc=rc)
            call LIS_verify(rc, &
            "Noah-MP.4.0.1 initial water in aquifer and saturated soil:"//&
                 " not defined")
        enddo

        ! leaf area index
        call ESMF_ConfigFindLabel(LIS_config, &
             "Noah-MP.4.0.1 initial leaf area index:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_lai, rc=rc)
            call LIS_verify(rc, &
                 "Noah-MP.4.0.1 initial leaf area index:"//&
                 " not defined")
        enddo

        ! optional gecros crop

        !!! no need to initialize Gecros crop model state variables
        !!! if the Noah-MP Gecros crop option is not turned on. SW
       if (NOAHMP401_struc(1)%crop_opt.eq.2) then
          call ESMF_ConfigFindLabel(LIS_config, &
               "Noah-MP.4.0.1 initial optional gecros crop:", rc = rc)
          do n=1,LIS_rc%nnest
              do i=1, 60
                  call ESMF_ConfigGetAttribute(LIS_config, NOAHMP401_struc(n)%init_gecros_state(i), rc=rc)
              end do
              call LIS_verify(rc, &
                   "Noah-MP.4.0.1 initial optional gecros crop:"//&
                   " not defined")
          enddo
        endif
    endif
     
    deallocate(nids)

 33   format(a47,i4)

    write(LIS_logunit, *) &
        "[INFO] Finish reading LIS configuration file for Noah-MP.4.0.1"
     
end subroutine NOAHMP401_readcrd
