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
! !ROUTINE: Crocus81_readcrd
! \label{Crocus81\_readcrd}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   10/18/19 : Mahdi Navari, Shugong Wang, initial implementation for LIS 7 and Crocus81
!   21 Jan 2021: Mahdi Navari, ground temperature removed form the lis.config, 
!                        for the stand-alone version, the value of TG was set to 273.15 in the Crocus81_main.F90
!   
! !INTERFACE:
subroutine Crocus81_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit , LIS_verify
    use Crocus81_lsmMod, only       : CROCUS81_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to Crocus81 model from
!  the LIS configuration file.
!
!EOP
    implicit none

    integer      :: rc 
    integer      :: n, i
    character*10 :: time 
    character*6  :: str_i
 
    write(LIS_logunit, *) "Start reading LIS configuration file for CROCUS81 model"
    
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "CROCUS81 model timestep: not defined")
        call LIS_parseTimeString(time, CROCUS81_struc(n)%ts)
    enddo
    
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 restart output interval:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc,"CROCUS81 restart output interval: not defined")
        call LIS_parseTimeString(time, CROCUS81_struc(n)%rstInterval)
    enddo
    
    !---------------------------!
    ! Constant Parameters       !
    !---------------------------!
    ! number of snow layer
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 nsnow:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%nsnow, rc=rc)
        call LIS_verify(rc, "CROCUS81 nsnow: not defined")
    enddo
 
    ! allocate memory for init_SNOWSWE using nsnow as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%init_SNOWSWE(CROCUS81_struc(n)%nsnow))
    enddo
    ! allocate memory for init_SNOWRHO using nsnow as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%init_SNOWRHO(CROCUS81_struc(n)%nsnow))
    enddo
    ! allocate memory for init_SNOWHEAT using nsnow as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%init_SNOWHEAT(CROCUS81_struc(n)%nsnow))
    enddo
    ! allocate memory for init_SNOWGRAN1 using nsnow as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%init_SNOWGRAN1(CROCUS81_struc(n)%nsnow))
    enddo
    ! allocate memory for init_SNOWGRAN2 using nsnow as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%init_SNOWGRAN2(CROCUS81_struc(n)%nsnow))
    enddo
    ! allocate memory for init_SNOWHIST using nsnow as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%init_SNOWHIST(CROCUS81_struc(n)%nsnow))
    enddo
    ! allocate memory for init_SNOWAGE using nsnow as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%init_SNOWAGE(CROCUS81_struc(n)%nsnow))
    enddo
    ! allocate memory for init_SNOWLIQ using nsnow as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%init_SNOWLIQ(CROCUS81_struc(n)%nsnow))
    enddo
    ! allocate memory for init_SNOWTEMP using nsnow as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%init_SNOWTEMP(CROCUS81_struc(n)%nsnow))
    enddo
    ! allocate memory for init_SNOWDZ using nsnow as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%init_SNOWDZ(CROCUS81_struc(n)%nsnow))
    enddo
 
    ! number of impurtites
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 nimpur:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%nimpur, rc=rc)
        call LIS_verify(rc, "CROCUS81 nimpur: not defined")
    enddo
 
    ! allocate memory for IMPWET using nimpur as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%IMPWET(CROCUS81_struc(n)%nimpur))
    enddo
    ! allocate memory for IMPDRY using nimpur as dimension
    do n=1, LIS_rc%nnest
        allocate(CROCUS81_struc(n)%IMPDRY(CROCUS81_struc(n)%nimpur))
    enddo
 
    ! SNOWRES_opt  = ISBA-SNOW3L turbulant exchange option
    !   'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
    !   'RIL' = Limit Richarson number under very stable
    ! conditions (currently testing)
    !   'M98'  = Martin et Lejeune 1998 : older computation for turbulent fluxes coefficents in Crocus
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWRES_opt:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWRES_opt, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWRES_opt: not defined")
    enddo
 
    ! True = coupled to MEB. surface fluxes are IMPOSED
    ! as an upper boundary condition to the explicit snow schemes. 
    ! If = False, then energy budget and fluxes are computed herein.
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 OMEB_BOOL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%OMEB_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 OMEB_BOOL: not defined")
    enddo
 

      ! True = Over permanent snow and ice, initialise WGI=WSAT, Hsnow>=10m and allow 0.8<SNOALB<0.85
      ! False = No specific treatment 
      call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 GLACIER_BOOL:", rc = rc) 
      do n=1, LIS_rc%nnest 
          call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%GLACIER_BOOL, rc=rc) 
          call LIS_verify(rc, "CROCUS81 GLACIER_BOOL: not defined") 
      enddo

    ! wind implicitation option  'OLD' = direct , 'NEW' = Taylor serie, order 1
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 HIMPLICIT_WIND_opt:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%HIMPLICIT_WIND_opt, rc=rc)
        call LIS_verify(rc, "CROCUS81 HIMPLICIT_WIND_opt: not defined")
    enddo
 
    ! time step of the integration
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 PTSTEP:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%PTSTEP, rc=rc)
        call LIS_verify(rc, "CROCUS81 PTSTEP: not defined")
    enddo

# if 0
    ! MN: For now we assume there is no energy transfer between snow and soil by setting surface 
    !     soil temperature to 273.15 in the lis.config in feature we will use surface soil temperature from an LSM.  
    !     Surface soil temperature (effective temperature the of layer lying below snow) (K)  (for snowcro.F90 
    !     we only use the surface layer ZP_TG(:,1))  (#nsoil depends on 2-L, 3-L DIF)

    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 TG:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%TG, rc=rc)
        call LIS_verify(rc, "CROCUS81 TG: not defined")
    enddo
# endif
    ! reference height of the wind
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 UREF:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%UREF, rc=rc)
        call LIS_verify(rc, "CROCUS81 UREF: not defined")
    enddo
 
    ! Reference height of the first atmospheric level (m)
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 ZREF:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%ZREF, rc=rc)
        call LIS_verify(rc, "CROCUS81 ZREF: not defined")
    enddo
 
    ! grid box average roughness length (m) (roughness length for momentum)
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 Z0NAT:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%Z0NAT, rc=rc)
        call LIS_verify(rc, "CROCUS81 Z0NAT: not defined")
    enddo
 
    ! roughness length for momentum (modd_diagn.F90 effective roughness length for heat(!?))
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 Z0EFF:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%Z0EFF, rc=rc)
        call LIS_verify(rc, "CROCUS81 Z0EFF: not defined")
    enddo
 
    ! grid box average roughness length for heat
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 Z0HNAT:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%Z0HNAT, rc=rc)
        call LIS_verify(rc, "CROCUS81 Z0HNAT: not defined")
    enddo


    ! if usemonalb == .true., then the alb value passed to lsmcrocus will be used as the background snow-free albedo term.  
    ! if usemonalb == .false., then alb will be set to 0.2 
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 use monthly albedo map:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%use_monthly_albedo_map, rc=rc)
        call LIS_verify(rc, "CROCUS81 use monthly albedo map: not defined")
    enddo


    ! Assumed first soil layer thickness (m)
    ! Used to calculate ground/snow heat flux   (D_G(:,1))
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 D_G:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%D_G, rc=rc)
        call LIS_verify(rc, "CROCUS81 D_G: not defined")
    enddo


    ! Mechanical transformation of snow grain and compaction + effect of wind on falling snow properties
    !	'NONE': No snowdrift scheme
    !	'DFLT': falling snow falls as purely dendritic
    ! 	'GA01': Gallee et al 2001
    !	'VI13': Vionnet et al 2013
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWDRIFT_opt:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWDRIFT_opt, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWDRIFT_opt: not defined")
    enddo
 
    ! Logicals for snowdrift sublimation
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWDRIFT_SUBLIM_BOOL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWDRIFT_SUBLIM_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWDRIFT_SUBLIM_BOOL: not defined")
    enddo
 
    ! Activate parametrization of solar absorption for polar regions (If True modify solar absorption as a function of solar zenithal angle)
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOW_ABS_ZENITH_BOOL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOW_ABS_ZENITH_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOW_ABS_ZENITH_BOOL: not defined")
    enddo
 
    ! Metamorphism scheme: B92 (historical version, Brun et al 92), C13, T07, F06 (see Carmagnola et al 2014)
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWMETAMO_opt:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWMETAMO_opt, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWMETAMO_opt: not defined")
    enddo
 
    ! Radiative transfer scheme. HSNOWRAD=B92 Brun et al 1992.  HSNOWRAD=T17 (Tuzet et al. 2017) (Libois et al. 2013) TARTES with impurities content scheme
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWRAD_opt:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWRAD_opt, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWRAD_opt: not defined")
    enddo
 
    ! Activate atmotartes scheme  (default=.FALSE. # This option is not stable yet, but it is supposed to compute the direct/diffuse ratio directly from atmospheric informations (AOD, Ozone column, Water column..))
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 ATMORAD_BOOL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%ATMORAD_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 ATMORAD_BOOL: not defined")
    enddo
 
    ! [init_surf_atmn.F90   --> wet deposit coefficient for each impurity type (g/m_/s)  ,   snowcro.F90 --> Dry and wet deposit coefficient from Forcing File(g/m_/s)   (# 1553 You can either feed the model with prescribed and constant deposition fluxes or introduce a wet and dry deposition field directly in the forcing file. )
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 IMPWET:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, CROCUS81_struc(n)%nimpur
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%IMPWET(i), rc=rc)
            call LIS_verify(rc, 'CROCUS81 IMPWET: not defined')
        enddo
    enddo
 
    ! [init_surf_atmn.F90   --> wet deposit coefficient for each impurity type (g)  ,   snowcro.F90 --> Dry and wet deposit coefficient from Forcing File(g/m_/s)
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 IMPDRY:", rc = rc)
    do n=1, LIS_rc%nnest
        do i = 1, CROCUS81_struc(n)%nimpur
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%IMPDRY(i), rc=rc)
            call LIS_verify(rc, 'CROCUS81 IMPDRY: not defined')
        enddo
    enddo
 
    ! New options for multiphysics version (Cluzet et al 2016). Falling snow scheme: V12 (Vionnet et al. 2012) , A76 (Anderson 1976), S02 (Lehning and al. 2002), P75 (Pahaut 1975)
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWFALL_opt:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWFALL_opt, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWFALL_opt: not defined")
    enddo
 
    ! Thermal conductivity scheme: Y81 (Yen 1981), I02 (Boone et al. 2002) C11 (Calonne et al. 2011)
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWCOND_opt:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWCOND_opt, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWCOND_opt: not defined")
    enddo
 
    ! liquid water content scheme: B92 (Brun et al. 1992) O04 (Oleson et al., 2004) S02 (SNOWPACK, Lehning et al, 2002) B02 (ISBA_ES, Boone et al. 2002)
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWHOLD_opt:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWHOLD_opt, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWHOLD_opt: not defined")
    enddo
 
    ! B92 snow compaction basis version and B93 for slightly different parameters   (NOTE: IN THE CODE S14, T11, )
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWCOMP_opt:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWCOMP_opt, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWCOMP_opt: not defined")
    enddo
 
    ! reference height is constant or variable from the snow surface: CST (constant from snow surface, i.e. Col de Porte) or VAR (variable from snow surface = snow depth has to be removed from reference height)
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWZREF_opt:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWZREF_opt, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWZREF_opt: not defined")
    enddo
 
    ! Snowmaking and Grooming options
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWCOMPACT_BOOL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWCOMPACT_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWCOMPACT_BOOL: not defined")
    enddo
 
    ! Snowmaking and Grooming options
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWMAK_BOOL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWMAK_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWMAK_BOOL: not defined")
    enddo
 
    ! Snowmaking and Grooming options
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWTILLER_BOOL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWTILLER_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWTILLER_BOOL: not defined")
    enddo
 
    ! Snowmaking and Grooming options
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SELF_PROD_BOOL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SELF_PROD_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 SELF_PROD_BOOL: not defined")
    enddo
 
    ! Snowmaking and Grooming options
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SNOWMAK_PROP_BOOL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SNOWMAK_PROP_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 SNOWMAK_PROP_BOOL: not defined")
    enddo
 
    ! Snowmaking and Grooming options
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 PRODSNOWMAK_BOOL:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%PRODSNOWMAK_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 PRODSNOWMAK_BOOL: not defined")
    enddo
 
    ! Use ASPECT (SLOPE DIRECTION) from LDT output    
    !  !Typical slope aspect in the grid  (clockwise from N)  
    !      call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 SLOPE_DIR:", rc = rc)  
    !      do n=1, LIS_rc%nnest 
    !          call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%SLOPE_DIR, rc=rc)
    !          call LIS_verify(rc, "CROCUS81 SLOPE_DIR: not defined") 
    !      enddo 

    ! Boolean option to partition total precipitation into snowfall and rainfall using Jordan 1991
    call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 boolean option to partition total precip:", rc = rc)
    do n=1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%Partition_total_precip_BOOL, rc=rc)
        call LIS_verify(rc, "CROCUS81 boolean option to partition total precip: not defined")
    enddo

    ! The following lines hard code the LDT NetCDF variable names. 
    do n=1, LIS_rc%nnest
        !CROCUS81_struc(n)%LDT_ncvar_GLACIER_BOOL = 'CROCUS81_GLACIER_BOOL'
        !CROCUS81_struc(n)%LDT_ncvar_TG       = 'Crocus_TG' !'CROCUS81_TG'
        CROCUS81_struc(n)%LDT_ncvar_SLOPE    = 'SLOPE' !'CROCUS81_SLOPE'
        CROCUS81_struc(n)%LDT_ncvar_ALB      = 'ALBEDO'
        !CROCUS81_struc(n)%LDT_ncvar_SOILCOND = 'Crocus_SOILCOND' !'CROCUS81_SOILCOND'
        CROCUS81_struc(n)%LDT_ncvar_PERMSNOWFRAC =  'GLACIERFRAC' !  'CROCUS81_PERMSNOWFRAC'
        CROCUS81_struc(n)%LDT_ncvar_SLOPE_DIR = 'ASPECT'   ! 'CROCUS81_SLOPE_DIR'
        CROCUS81_struc(n)%LDT_ncvar_SAND     = 'SAND' ! 'CROCUS81_SAND'
        CROCUS81_struc(n)%LDT_ncvar_SILT     = 'SILT'
        CROCUS81_struc(n)%LDT_ncvar_CLAY     = 'CLAY'
        CROCUS81_struc(n)%LDT_ncvar_POROSITY = 'POROSITY'
    enddo

    ! set default restart format to netcdf
    do n=1,LIS_rc%nnest
        CROCUS81_struc(n)%rformat = "netcdf"
    enddo
    ! restart run, read restart file
    if (trim(LIS_rc%startcode) == "restart") then 
        Call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 restart file:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%rfile, rc=rc)
            call LIS_verify(rc, "CROCUS81 restart file: not defined")
        enddo
        
        Call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 restart file format:", rc=rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%rformat, rc=rc)
            call LIS_verify(rc, "CROCUS81 restart file format: not defined")
        enddo
    ! cold start run, read initial state variables
    else 
        ! Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWSWE:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, CROCUS81_struc(n)%nsnow
                call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWSWE(i), rc=rc)
            end do
            call LIS_verify(rc, "CROCUS81 initial SNOWSWE: not defined")
        enddo

        ! Snow layer(s) averaged density (kg/m3)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWRHO:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, CROCUS81_struc(n)%nsnow
                call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWRHO(i), rc=rc)
            end do
            call LIS_verify(rc, "CROCUS81 initial SNOWRHO: not defined")
        enddo

        ! Snow layer(s) heat content (J/m2)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWHEAT:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, CROCUS81_struc(n)%nsnow
                call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWHEAT(i), rc=rc)
            end do
            call LIS_verify(rc, "CROCUS81 initial SNOWHEAT: not defined")
        enddo

        ! snow surface albedo
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWALB:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWALB, rc=rc)
            call LIS_verify(rc, "CROCUS81 initial SNOWALB: not defined")
        enddo

        ! Snow layers grain feature 1
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWGRAN1:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, CROCUS81_struc(n)%nsnow
                call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWGRAN1(i), rc=rc)
            end do
            call LIS_verify(rc, "CROCUS81 initial SNOWGRAN1: not defined")
        enddo

        ! Snow layer grain feature 2
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWGRAN2:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, CROCUS81_struc(n)%nsnow
                call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWGRAN2(i), rc=rc)
            end do
            call LIS_verify(rc, "CROCUS81 initial SNOWGRAN2: not defined")
        enddo

        ! Snow layer grain historical parameter (only for non dendritic snow) (-) in {0-5}
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWHIST:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, CROCUS81_struc(n)%nsnow
                call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWHIST(i), rc=rc)
            end do
            call LIS_verify(rc, "CROCUS81 initial SNOWHIST: not defined")
        enddo

        ! Age since snowfall (day)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWAGE:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, CROCUS81_struc(n)%nsnow
                call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWAGE(i), rc=rc)
            end do
            call LIS_verify(rc, "CROCUS81 initial SNOWAGE: not defined")
        enddo

        ! Snow layer(s) liquid water content (m)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWLIQ:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, CROCUS81_struc(n)%nsnow
                call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWLIQ(i), rc=rc)
            end do
            call LIS_verify(rc, "CROCUS81 initial SNOWLIQ: not defined")
        enddo

        ! Snow layer(s) temperature (K)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWTEMP:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, CROCUS81_struc(n)%nsnow
                call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWTEMP(i), rc=rc)
            end do
            call LIS_verify(rc, "CROCUS81 initial SNOWTEMP: not defined")
        enddo

        ! Snow layer(s) thickness (m)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWDZ:", rc = rc)
        do n=1,LIS_rc%nnest
            do i=1, CROCUS81_struc(n)%nsnow
                call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWDZ(i), rc=rc)
            end do
            call LIS_verify(rc, "CROCUS81 initial SNOWDZ: not defined")
        enddo

        ! Soil/snow interface heat flux (W/m2)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial GRNDFLUX:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_GRNDFLUX, rc=rc)
            call LIS_verify(rc, "CROCUS81 initial GRNDFLUX: not defined")
        enddo

        ! Blowing snow sublimation (kg/m2/s) NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNDRIFT:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNDRIFT, rc=rc)
            call LIS_verify(rc, "CROCUS81 initial SNDRIFT: not defined")
        enddo

        ! Richardson number (-)  NOTE: RI has not been initialized in CALL_MODEL (If not OMED initalized to undefined in the snow3L_isba.F90)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial RI_n:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_RI_n, rc=rc)
            call LIS_verify(rc, "CROCUS81 initial RI_n: not defined")
        enddo

        ! Drag coefficient for momentum over snow (-)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial CDSNOW:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_CDSNOW, rc=rc)
            call LIS_verify(rc, "CROCUS81 initial CDSNOW: not defined")
        enddo

        ! Friction velocity over snow (m/s);
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial USTARSNOW:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_USTARSNOW, rc=rc)
            call LIS_verify(rc, "CROCUS81 initial USTARSNOW: not defined")
        enddo

        ! Drag coefficient for heat over snow  (-)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial CHSNOW:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_CHSNOW, rc=rc)
            call LIS_verify(rc, "CROCUS81 initial CHSNOW: not defined")
        enddo

        ! Snowmaking thickness (m)
        call ESMF_ConfigFindLabel(LIS_config, "CROCUS81 initial SNOWMAK_dz:", rc = rc)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, CROCUS81_struc(n)%init_SNOWMAK_dz, rc=rc)
            call LIS_verify(rc, "CROCUS81 initial SNOWMAK_dz: not defined")
        enddo

    end if
     
    write(LIS_logunit, *) "Finish reading LIS configuration file for CROCUS81 model"
     
end subroutine CROCUS81_readcrd
