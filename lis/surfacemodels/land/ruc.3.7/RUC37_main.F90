!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: RUC37_main
! \label{RUC37_main}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   1/15/15: Shugong Wang; initial implementation for RUC37 with LIS-7
!
! !INTERFACE:
subroutine RUC37_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod 
    use RUC37_lsmMod
   !use other modules
    use ruc_debug 
    implicit none
! !ARGUMENTS:
    integer, intent(in)  :: n
    integer              :: t
    integer              :: i
    real                 :: dt
    real                 :: lat, lon
    integer              :: row, col
    integer              :: year, month, day, hour, minute, second
    logical              :: alarmCheck

!
! !DESCRIPTION:
!  This is the entry point for calling the RUC37 physics.
!  This routine calls the {\tt ruc_driver_37 } routine that performs the
!  land surface computations, to solve for water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for RUC37
    integer              :: tmp_year               ! year [-]
    integer              :: tmp_month              ! month [-]
    integer              :: tmp_day                ! day [-]
    integer              :: tmp_hour               ! hour [-]
    integer              :: tmp_minute             ! minute [-]
    integer              :: tmp_ktime              ! time step 
    real                 :: tmp_lwdown             ! downward longwave radiation flux at surface (w m-2) [forcing] [W m-2]
    real                 :: tmp_swdown             ! downward shortwave radiation flux at surface (w m-2) [forcing] [W m-2]
    real                 :: tmp_psurf              ! surface atmospheric pressure (pa) [forcing] [Pa]
    real                 :: tmp_rainf              ! rainfall rate (kg m-2 s-1) [forcing] [kg m-2 s-1]
    real                 :: tmp_snowf              ! snowfall rate (kg m-2 s-1) [forcing] [Kg/m2s]
    real                 :: tmp_tair               ! air temperature (k) [forcing] [K]
    real                 :: tmp_qair               ! surface specific humidity (kg kg-1) [forcing] [kg kg-1]
    real                 :: tmp_wind_e             ! eastward wind speed (m s-1) [forcing] [m s-1]
    real                 :: tmp_wind_n             ! northward wind speed (m s-1) [forcing] [m s-1]
    integer              :: tmp_vegetype           ! vegetation category [-]
    integer              :: tmp_soiltype           ! soil category [-]
    real                 :: tmp_dt                 ! time step (seconds). [s]
    real, allocatable    :: tmp_soil_layer_thickness(:) ! thicknesses of each soil level (m) [m]
    real                 :: tmp_zlvl 
    real                 :: tmp_zlvl_wind 
    logical              :: tmp_use_local_param    ! .true. to use table values for albbck, shdfac, and z0brd; .false. to use values for albbck, shdfac, and z0brd as set in this driver routine [-]
    logical              :: tmp_use_2d_lai_map     ! if rdlai2d == .true., then the xlai value that we pass to lsmruc will be used. if rdlai2d == .false., then xlai will be computed within lsmruc, from table minimum and maximum values in vegparm.tbl, and the current green vegetation fraction. [-]
    logical              :: tmp_use_monthly_albedo_map ! if usemonalb == .true., then the alb value passed to lsmruc will be used as the background snow-free albedo term.  if usemonalb == .false., then alb will be computed within lsmruc from minimum and maximum values in vegparm.tbl, and the current green vegetation fraction. [-]
    integer              :: tmp_option_iz0tlnd     ! option to turn on (iz0tlnd=1) or off (iz0tlnd=0) the vegetation-category-dependent calculation of the zilitinkivich coefficient czil in the sfcdif subroutines. [-]
    integer              :: tmp_option_sfcdif      ! option to use previous (sfcdif_option=0) or updated (sfcdif_option=1) version of sfcdif subroutine. [-]
    character(len=256)   :: tmp_landuse_tbl_name   ! noah model landuse parameter table [-]
    character(len=256)   :: tmp_soil_tbl_name      ! noah model soil parameter table [-]
    character(len=256)   :: tmp_gen_tbl_name       ! noah model soil parameter table [-]
    character(len=256)   :: tmp_landuse_scheme_name ! landuse classification scheme [-]
    character(len=256)   :: tmp_soil_scheme_name   ! soil classification scheme [-]
    integer              :: tmp_nsoil              ! number of soil levels. [-]
    integer              :: tmp_water_class_num    ! number of water category in llanduse classification [-]
    integer              :: tmp_ice_class_num      ! number of ice category in llanduse classification [-]
    integer              :: tmp_urban_class_num    ! number of urban category in llanduse classification [-]
    real, allocatable    :: tmp_albedo_monthly(:)  ! monthly values of background (i.e., snow-free) albedo ( fraction [0.0-1.0] ) [-]
    real, allocatable    :: tmp_shdfac_monthly(:)  ! monthly values for green vegetation fraction ( fraction [0.0-1.0] ) [-]
    real, allocatable    :: tmp_z0brd_monthly(:)   ! monthly values for background (i.e., snow-free) roughness length ( m ) [m]
    real, allocatable    :: tmp_lai_monthly(:)     ! monthly values for leaf area index ( dimensionless ) [-]
    real                 :: tmp_albbck             ! background snow-free albedo (0.0-1.0). [-]
    real                 :: tmp_tbot               ! deep-soil time-invariant temperature (k).  representing sort of a mean annual air temperature. [K]
    real                 :: tmp_snoalb             ! maximum snow albedo over deep snow (0.0-1.0) [-]
    real                 :: tmp_emiss              ! surface emissivity (0.0 - 1.0). [-]
    real                 :: tmp_ch                 ! exchange coefficient for head and moisture (m s-1). [s/m]
    real                 :: tmp_cm                 ! exchange coefficient for momentum (m s-1). [s/m]
    real                 :: tmp_sneqv              ! water equivalent of accumulated snow depth (m). [m]
    real                 :: tmp_snowh              ! physical snow depth (m). [m]
    real                 :: tmp_snowc              ! fractional snow cover ( fraction [0.0-1.0] ) [-]
    real                 :: tmp_canwat             ! canopy moisture content (kg m-2) [kg m-2]
    real                 :: tmp_alb                ! surface albedo including possible snow-cover effect.  this is set in lsmruc, [-]
    real, allocatable    :: tmp_smc(:)             ! total soil moisture content (m3 m-3) [m^3 m-3]
    real, allocatable    :: tmp_sho(:)             ! liquid soil moisture content (m3 m-3) [m^3 m-3]
    real, allocatable    :: tmp_stc(:)             ! soil temperature (k) [K]
    real, allocatable    :: tmp_smfr(:)            ! soil ice content (m3 m-3) [m^3 m-3]
    real, allocatable    :: tmp_keepfr(:)          ! soil ice flag: 0 or 1
    real                 :: tmp_tskin              ! skin temperature (k) [K]
    real                 :: tmp_qvg                ! mixing ratio at the surface ( kg kg{-1} ) [kg kg-1]
    real                 :: tmp_qsfc               ! specific humidity at the surface ( kg kg{-1} ) [kg kg-1]
    real                 :: tmp_qcg                ! cloud water mixing ratio at the surface ( kg kg{-1} ) [kg kg-1]
    real                 :: tmp_qsg                ! surface water vapor mixing ratio at satration (kg kg-1) [kg/kg]
    real                 :: tmp_snt75cm            ! snow temperature at 7.5 cm depth (k) [K]
    real                 :: tmp_tsnav              ! average snow temperature in k [K]
    real                 :: tmp_soilm              ! total soil column moisture content, frozen and unfrozen ( m ) [m]
    real                 :: tmp_smroot             ! available soil moisture in the root zone ( fraction [smcwlt-smcmax] [m^3 m-3]
    real                 :: tmp_shdfac             ! shading factor ( - )
    real                 :: tmp_shdmin             ! min shading factor ( - )
    real                 :: tmp_shdmax             ! max shading factor ( - )
    real                 :: tmp_qh                 ! sensible heat flux ( w m{-2} ) [W m-2]
    real                 :: tmp_qle                ! latent heat flux (evapotranspiration) ( w m{-2} ) [W m-2]
    real                 :: tmp_eta                ! latent heat flux (evapotranspiration) ( kg m{-2} s{-1} ) [kg m-2 s-1]
    real                 :: tmp_qg                 ! soil heat flux ( w m{-2} ) [W m-2]
    real                 :: tmp_rhosnf             ! density of frozen precipitation (kg m{-3}) [kg/m3]
    real                 :: tmp_precipfr           ! time-step frozen precipitation (kg m{-2}) [kg m-2]
    real                 :: tmp_snowfallac         ! run total snowfall accumulation (kg m{-2}) [kg m-2]
    real                 :: tmp_snthresh           ! snow depth threshold [m]
    real                 :: tmp_acsnow             ! run total frozen precipitation (kg m{-2}) [kg m-2]
    real                 :: tmp_sfcevp             ! run total evaporation flux  (kg m{-2}) [kg m-2]
    real                 :: tmp_snomlt             ! snow melt water ( m ) [m]
    real                 :: tmp_dew                ! dewfall (or frostfall for t<273.15) ( m ) [m]
    real                 :: tmp_drip               ! throughfall of precipitation from canopy (kg m{-2} s{-1}) [kg m-2 s-1]
    real                 :: tmp_qs                 ! surface runoff, not infiltrating the soil ( m s{-1} ) [kg m-2 s-1]
    real                 :: tmp_qsb                ! subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} ) [kg m-2 s-1]
    real                 :: tmp_snflx              ! snow heat flux (w/m^2: negative, if downward from surface) [W m-2]
    real                 :: tmp_edir               ! latent heat flux component: direct soil evaporation ( w m{-2} ) [W m-2]
    real                 :: tmp_ec                 ! latent heat flux component: canopy water evaporation ( w m{-2} ) [W m-2]
    real                 :: tmp_ett                ! latent heat flux component: total plant transpiration ( w m{-2} ) [W m-2]
    real                 :: tmp_esnow              ! sublimation from snowpack (w m{-2}) [W m-2]
    real                 :: tmp_qmax
    real                 :: tmp_qmin
    real                 :: tmp_psis
    real                 :: tmp_ksat
    real                 :: tmp_bclh
    real                 :: tmp_qwrtz
    real                 :: tmp_wilt
    real                 :: tmp_ref

    allocate( tmp_soil_layer_thickness( RUC37_struc(n)%nsoil ) )
    allocate( tmp_albedo_monthly( 12 ) )
    allocate( tmp_shdfac_monthly( 12 ) )
    allocate( tmp_z0brd_monthly( 12 ) )
    allocate( tmp_lai_monthly( 12 ) )
    allocate( tmp_smc   ( RUC37_struc(n)%nsoil ) )
    allocate( tmp_sho   ( RUC37_struc(n)%nsoil ) )
    allocate( tmp_stc   ( RUC37_struc(n)%nsoil ) )
    allocate( tmp_smfr  ( RUC37_struc(n)%nsoil ) )
    allocate( tmp_keepfr( RUC37_struc(n)%nsoil ) )

    ! check RUC37 alarm. If alarm is ring, run model. 
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "RUC37 model alarm")
    if (alarmCheck) Then
        do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            dt = LIS_rc%ts
            row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
            lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
            debug_row = row
            debug_col = col 
            tmp_ktime = LIS_rc%tscount(n) 

            ! retrieve forcing data from RUC37_struc(n)%ruc37(t) and assign to local variables
            ! lwdown: downward longwave radiation flux at surface (w m-2) [forcing]
            tmp_lwdown     = RUC37_struc(n)%ruc37(t)%lwdown / RUC37_struc(n)%forc_count
 
            ! swdown: downward shortwave radiation flux at surface (w m-2) [forcing]
            tmp_swdown     = RUC37_struc(n)%ruc37(t)%swdown / RUC37_struc(n)%forc_count
 
            ! psurf: surface atmospheric pressure (pa) [forcing]
            tmp_psurf      = RUC37_struc(n)%ruc37(t)%psurf  / RUC37_struc(n)%forc_count
 
            ! rainf: rainfall rate (kg m-2 s-1) [forcing]
            tmp_rainf      = RUC37_struc(n)%ruc37(t)%rainf  / RUC37_struc(n)%forc_count
 
            ! snowf: snowfall rate (kg m-2 s-1) [forcing]
            if(LIS_Forc_snowf%selectOpt .eq. 1) then 
                tmp_snowf  = RUC37_struc(n)%ruc37(t)%snowf  / RUC37_struc(n)%forc_count
            endif
 
            ! tair: air temperature (k) [forcing]
            tmp_tair       = RUC37_struc(n)%ruc37(t)%tair   / RUC37_struc(n)%forc_count
 
            ! qair: surface specific humidity (kg kg-1) [forcing]
            tmp_qair       = RUC37_struc(n)%ruc37(t)%qair   / RUC37_struc(n)%forc_count
 
            ! wind_e: eastward wind speed (m s-1) [forcing]
            tmp_wind_e     = RUC37_struc(n)%ruc37(t)%wind_e / RUC37_struc(n)%forc_count
 
            ! wind_n: northward wind speed (m s-1) [forcing]
            tmp_wind_n     = RUC37_struc(n)%ruc37(t)%wind_n / RUC37_struc(n)%forc_count
 
            ! 
            ! check validity of lwdown
            if(tmp_lwdown .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable lwdown in RUC37"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of swdown
            if(tmp_swdown .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable swdown in RUC37"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of psurf
            if(tmp_psurf .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable psurf in RUC37"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of rainf
            if(tmp_rainf .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable rainf in RUC37"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of snowf
            if((LIS_Forc_snowf%selectOpt .eq. 1) .and. (tmp_snowf .eq. LIS_rc%udef)) then
                write(LIS_logunit, *) "undefined value found for forcing variable snowf in RUC37"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of tair
            if(tmp_tair .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable tair in RUC37"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of qair
            if(tmp_qair .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable qair in RUC37"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of wind_e
            if(tmp_wind_e .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable wind_e in RUC37"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of wind_n
            if(tmp_wind_n .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable wind_n in RUC37"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! 
            tmp_year   = LIS_rc%yr
            tmp_month  = LIS_rc%mo
            tmp_day    = LIS_rc%da
            tmp_hour   = LIS_rc%hr
            tmp_minute = LIS_rc%mn
 
            ! get parameters 
            tmp_vegetype                            = RUC37_struc(n)%ruc37(t)%vegetype                        
            tmp_soiltype                            = RUC37_struc(n)%ruc37(t)%soiltype                        
            tmp_dt                                  = RUC37_struc(n)%dt                              
            tmp_soil_layer_thickness(:)             = RUC37_struc(n)%soil_layer_thickness(:)            
            tmp_zlvl                                = RUC37_struc(n)%zlvl
            tmp_zlvl_wind                           = RUC37_struc(n)%zlvl_wind
            tmp_use_local_param                     = RUC37_struc(n)%use_local_param                 
            tmp_use_2d_lai_map                      = RUC37_struc(n)%use_2d_lai_map                  
            tmp_use_monthly_albedo_map              = RUC37_struc(n)%use_monthly_albedo_map          
            tmp_option_iz0tlnd                      = RUC37_struc(n)%option_iz0tlnd                  
            tmp_option_sfcdif                       = RUC37_struc(n)%option_sfcdif                   
            tmp_landuse_tbl_name                    = RUC37_struc(n)%landuse_tbl_name                
            tmp_soil_tbl_name                       = RUC37_struc(n)%soil_tbl_name                   
            tmp_gen_tbl_name                        = RUC37_struc(n)%gen_tbl_name                   
            tmp_landuse_scheme_name                 = RUC37_struc(n)%landuse_scheme_name             
            tmp_soil_scheme_name                    = RUC37_struc(n)%soil_scheme_name                
            tmp_nsoil                               = RUC37_struc(n)%nsoil                           
            tmp_water_class_num                     = RUC37_struc(n)%water_class_num                 
            tmp_ice_class_num                       = RUC37_struc(n)%ice_class_num                   
            tmp_urban_class_num                     = RUC37_struc(n)%urban_class_num                   
            tmp_albedo_monthly(:)                   = RUC37_struc(n)%ruc37(t)%albedo_monthly(:)                  
            tmp_shdfac_monthly(:)                   = RUC37_struc(n)%ruc37(t)%shdfac_monthly(:)                  
!            tmp_z0brd_monthly(:)                    = RUC37_struc(n)%ruc37(t)%z0brd_monthly(:)                   
            tmp_lai_monthly(:)                      = RUC37_struc(n)%ruc37(t)%lai_monthly(:)                     
            tmp_albbck                              = RUC37_struc(n)%ruc37(t)%albbck                          
            tmp_tbot                                = RUC37_struc(n)%ruc37(t)%tbot                            
            tmp_snoalb                              = RUC37_struc(n)%ruc37(t)%snoalb                          
 
            ! get state variables
            tmp_emiss      = RUC37_struc(n)%ruc37(t)%emiss  
            tmp_ch         = RUC37_struc(n)%ruc37(t)%ch     
            tmp_cm         = RUC37_struc(n)%ruc37(t)%cm     
            tmp_sneqv      = RUC37_struc(n)%ruc37(t)%sneqv  
            tmp_snowh      = RUC37_struc(n)%ruc37(t)%snowh  
            tmp_snowc      = RUC37_struc(n)%ruc37(t)%snowc  
            tmp_canwat     = RUC37_struc(n)%ruc37(t)%canwat 
            tmp_alb        = RUC37_struc(n)%ruc37(t)%alb    
            tmp_smc(:)     = RUC37_struc(n)%ruc37(t)%smc(:)    
            tmp_sho(:)     = RUC37_struc(n)%ruc37(t)%sho(:)    
            tmp_stc(:)     = RUC37_struc(n)%ruc37(t)%stc(:)    
            tmp_smfr(:)    = RUC37_struc(n)%ruc37(t)%smfr(:)    
            tmp_keepfr(:)  = RUC37_struc(n)%ruc37(t)%keepfr(:)    
            tmp_tskin      = RUC37_struc(n)%ruc37(t)%tskin  
            tmp_qvg        = RUC37_struc(n)%ruc37(t)%qvg    
            tmp_qsfc       = RUC37_struc(n)%ruc37(t)%qsfc    
            tmp_qcg        = RUC37_struc(n)%ruc37(t)%qcg    
            tmp_qsg        = RUC37_struc(n)%ruc37(t)%qsg    
            tmp_snt75cm    = RUC37_struc(n)%ruc37(t)%snt75cm
            tmp_tsnav      = RUC37_struc(n)%ruc37(t)%tsnav 
            tmp_soilm      = RUC37_struc(n)%ruc37(t)%soilm  
            tmp_smroot     = RUC37_struc(n)%ruc37(t)%smroot 
 
            tmp_snowfallac   =  RUC37_struc(n)%ruc37(t)%snowfallac  
            tmp_snthresh     =  RUC37_struc(n)%ruc37(t)%snthresh
            tmp_acsnow       =  RUC37_struc(n)%ruc37(t)%acsnow      
            tmp_sfcevp       =  RUC37_struc(n)%ruc37(t)%sfcevp      

            ! call model physics 
            call ruc_driver_37(tmp_year              , & ! in    - year [-]
                               tmp_month             , & ! in    - month [-]
                               tmp_day               , & ! in    - day [-]
                               tmp_hour              , & ! in    - hour [-]
                               tmp_minute            , & ! in    - minute [-]
                               tmp_ktime             , & ! in    - time step [-] 
                               tmp_lwdown            , & ! in    - downward longwave radiation flux at surface (w m-2) [forcing] [W m-2]
                               tmp_swdown            , & ! in    - downward shortwave radiation flux at surface (w m-2) [forcing] [W m-2]
                               tmp_psurf             , & ! in    - surface atmospheric pressure (pa) [forcing] [Pa]
                               tmp_rainf             , & ! in    - rainfall rate (kg m-2 s-1) [forcing] [kg m-2 s-1]
                               tmp_snowf             , & ! in    - snowfall rate (kg m-2 s-1) [forcing] [Kg/m2s]
                               tmp_tair              , & ! in    - air temperature (k) [forcing] [K]
                               tmp_qair              , & ! in    - surface specific humidity (kg kg-1) [forcing] [kg kg-1]
                               tmp_wind_e            , & ! in    - eastward wind speed (m s-1) [forcing] [m s-1]
                               tmp_wind_n            , & ! in    - northward wind speed (m s-1) [forcing] [m s-1]
                               tmp_vegetype          , & ! in    - vegetation category [-]
                               tmp_soiltype          , & ! in    - soil category [-]
                               tmp_dt                , & ! in    - time step (seconds). [s]
                               tmp_soil_layer_thickness, & ! in    - thicknesses of each soil level (m) [m]
                               tmp_zlvl              , &
                               tmp_zlvl_wind         , &
                               tmp_use_local_param   , & ! in    - .true. to use table values for albbck, shdfac, and z0brd; .false. to use values for albbck, shdfac, and z0brd as set in this driver routine [-]
                               tmp_use_2d_lai_map    , & ! in    - if rdlai2d == .true., then the xlai value that we pass to lsmruc will be used. if rdlai2d == .false., then xlai will be computed within lsmruc, from table minimum and maximum values in vegparm.tbl, and the current green vegetation fraction. [-]
                               tmp_use_monthly_albedo_map, & ! in    - if usemonalb == .true., then the alb value passed to lsmruc will be used as the background snow-free albedo term.  if usemonalb == .false., then alb will be computed within lsmruc from minimum and maximum values in vegparm.tbl, and the current green vegetation fraction. [-]
                               tmp_option_iz0tlnd    , & ! in    - option to turn on (iz0tlnd=1) or off (iz0tlnd=0) the vegetation-category-dependent calculation of the zilitinkivich coefficient czil in the sfcdif subroutines. [-]
                               tmp_option_sfcdif     , & ! in    - option to use previous (sfcdif_option=0) or updated (sfcdif_option=1) version of sfcdif subroutine. [-]
                               tmp_landuse_tbl_name  , & ! in    - noah model landuse parameter table [-]
                               tmp_soil_tbl_name     , & ! in    - noah model soil parameter table [-]
                               tmp_gen_tbl_name      , & ! in    - noah model general  parameter table [-]
                               tmp_landuse_scheme_name, & ! in    - landuse classification scheme [-]
                               tmp_soil_scheme_name  , & ! in    - soil classification scheme [-]
                               tmp_nsoil             , & ! in    - number of soil levels. [-]
                               tmp_water_class_num   , & ! in    - number of water category in llanduse classification [-]
                               tmp_ice_class_num     , & ! in    - number of ice category in llanduse classification [-]
                               tmp_urban_class_num   , & ! in    - number of urban category in llanduse classification [-]
                               tmp_albedo_monthly    , & ! in    - monthly values of background (i.e., snow-free) albedo ( fraction [0.0-1.0] ) [-]
                               tmp_shdfac_monthly    , & ! in    - monthly values for green vegetation fraction ( fraction [0.0-1.0] ) [-]
                               tmp_z0brd_monthly     , & ! in    - monthly values for background (i.e., snow-free) roughness length ( m ) [m]
                               tmp_lai_monthly       , & ! in    - monthly values for leaf area index ( dimensionless ) [-]
                               tmp_albbck            , & ! in    - background snow-free albedo (0.0-1.0). [-]
                               tmp_tbot              , & ! in    - deep-soil time-invariant temperature (k).  representing sort of a mean annual air temperature. [K]
                               tmp_snoalb            , & ! in    - maximum snow albedo over deep snow (0.0-1.0) [-]
                               tmp_emiss             , & ! inout - surface emissivity (0.0 - 1.0). [-]
                               tmp_ch                , & ! inout - exchange coefficient for head and moisture (m s-1). [s/m]
                               tmp_cm                , & ! inout - exchange coefficient for momentum (m s-1). [s/m]
                               tmp_sneqv             , & ! inout - water equivalent of accumulated snow depth (m). [m]
                               tmp_snowh             , & ! inout - physical snow depth (m). [m]
                               tmp_snowc             , & ! inout - fractional snow cover ( fraction [0.0-1.0] ) [-]
                               tmp_canwat            , & ! inout - canopy moisture content (kg m-2) [kg m-2]
                               tmp_alb               , & ! inout - surface albedo including possible snow-cover effect.  this is set in lsmruc, [-]
                               tmp_smc               , & ! inout - total soil moisture content (m3 m-3) [m^3 m-3]
                               tmp_sho               , & ! inout - liquid soil moisture content (m3 m-3) [m^3 m-3]
                               tmp_stc               , & ! inout - soil temperature (k) [K]
                               tmp_smfr              , & ! inout - soil ice content [m^3 m-3]
                               tmp_keepfr            , & ! inout - frozen soil flag
                               tmp_tskin             , & ! inout - skin temperature (k) [K]
                               tmp_qvg               , & ! inout - mixing ratio at the surface ( kg kg{-1} ) [kg kg-1]
                               tmp_qsfc              , & ! inout - specific humidity at the surface ( kg kg{-1} ) [kg kg-1]
                               tmp_qcg               , & ! inout - cloud water mixing ratio at the surface ( kg kg{-1} ) [kg kg-1]
                               tmp_qsg               , & ! inout - surface water vapor mixing ratio at satration (kg kg-1) [kg/kg]
                               tmp_snt75cm           , & ! inout - snow temperature at 7.5 cm depth (k) [K]
                               tmp_tsnav             , & ! inout - average snow temperature in k [K]
                               tmp_soilm             , & ! inout - total soil column moisture content, frozen and unfrozen ( m ) [m]
                               tmp_smroot            , & ! inout - available soil moisture in the root zone ( fraction [smcwlt-smcmax] [m^3 m-3]
                               tmp_shdfac            , & ! out   - shading factor ( 0. - 1. )
                               tmp_shdmin            , & ! out   - min shading factor ( 0. - 1. )
                               tmp_shdmax            , & ! out   - max shading factor ( 0. - 1. )
                               tmp_qh                , & ! out   - sensible heat flux ( w m{-2} ) [W m-2]
                               tmp_qle               , & ! out   - latent heat flux (evapotranspiration) ( w m{-2} ) [W m-2]
                               tmp_eta               , & ! out   - latent heat flux (evapotranspiration) ( kg m{-2} s{-1} ) [kg m-2 s-1]
                               tmp_qg                , & ! out   - soil heat flux ( w m{-2} ) [W m-2]
                               tmp_rhosnf            , & ! out   - density of frozen precipitation (kg m{-3}) [kg/m3]
                               tmp_precipfr          , & ! out   - time-step frozen precipitation (kg m{-2}) [kg m-2]
                               tmp_snowfallac        , & ! out   - run total snowfall accumulation (kg m{-2}) [kg m-2]
                               tmp_qmax              , & ! out
                               tmp_qmin              , & ! out
                               tmp_psis              , & ! out
                               tmp_ksat              , & ! out
                               tmp_bclh              , & ! out
                               tmp_qwrtz             , & ! out
                               tmp_wilt              , & ! out
                               tmp_ref               , & ! out
                               tmp_snthresh          , & ! out   - snow depth threshold [m]
                               tmp_acsnow            , & ! out   - run total frozen precipitation (kg m{-2}) [kg m-2]
                               tmp_sfcevp            , & ! out   - run total evaporation flux  (kg m{-2}) [kg m-2]
                               tmp_snomlt            , & ! out   - snow melt water ( m ) [m]
                               tmp_dew               , & ! out   - dewfall (or frostfall for t<273.15) ( m ) [m]
                               tmp_drip              , & ! out   - throughfall of precipitation from canopy (kg m{-2} s{-1}) [kg m-2 s-1]
                               tmp_qs                , & ! out   - surface runoff, not infiltrating the soil ( m s{-1} ) [kg m-2 s-1]
                               tmp_qsb               , & ! out   - subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} ) [kg m-2 s-1]
                               tmp_snflx             , & ! out   - snow heat flux (w/m^2: negative, if downward from surface) [W m-2]
                               tmp_edir              , & ! out   - latent heat flux component: direct soil evaporation ( w m{-2} ) [W m-2]
                               tmp_ec                , & ! out   - latent heat flux component: canopy water evaporation ( w m{-2} ) [W m-2]
                               tmp_ett               , & ! out   - latent heat flux component: total plant transpiration ( w m{-2} ) [W m-2]
                               tmp_esnow             )   ! out   - sublimation from snowpack (w m{-2}) [W m-2]
 
    
            ! save state variables from local variables to global variables
            RUC37_struc(n)%ruc37(t)%emiss      = tmp_emiss     
            RUC37_struc(n)%ruc37(t)%ch         = tmp_ch        
            RUC37_struc(n)%ruc37(t)%cm         = tmp_cm        
            RUC37_struc(n)%ruc37(t)%sneqv      = tmp_sneqv     
            RUC37_struc(n)%ruc37(t)%snowh      = tmp_snowh     
            RUC37_struc(n)%ruc37(t)%snowc      = tmp_snowc     
            RUC37_struc(n)%ruc37(t)%canwat     = tmp_canwat    
            RUC37_struc(n)%ruc37(t)%alb        = tmp_alb       
            RUC37_struc(n)%ruc37(t)%albbck     = tmp_albbck       
            RUC37_struc(n)%ruc37(t)%smc(:)     = tmp_smc(:)    
            RUC37_struc(n)%ruc37(t)%sho(:)     = tmp_sho(:)    
            RUC37_struc(n)%ruc37(t)%stc(:)     = tmp_stc(:)    
            RUC37_struc(n)%ruc37(t)%smfr(:)    = tmp_smfr(:)    
            RUC37_struc(n)%ruc37(t)%keepfr(:)  = tmp_keepfr(:)    
            RUC37_struc(n)%ruc37(t)%tskin      = tmp_tskin     
            RUC37_struc(n)%ruc37(t)%qvg        = tmp_qvg       
            RUC37_struc(n)%ruc37(t)%qsfc       = tmp_qsfc       
            RUC37_struc(n)%ruc37(t)%qcg        = tmp_qcg       
            RUC37_struc(n)%ruc37(t)%qsg        = tmp_qsg       
            RUC37_struc(n)%ruc37(t)%snt75cm    = tmp_snt75cm   
            RUC37_struc(n)%ruc37(t)%tsnav      = tmp_tsnav
            RUC37_struc(n)%ruc37(t)%soilm      = tmp_soilm     
            RUC37_struc(n)%ruc37(t)%smroot     = tmp_smroot    
            RUC37_struc(n)%ruc37(t)%snowfallac    = tmp_snowfallac   
            RUC37_struc(n)%ruc37(t)%snthresh      = tmp_snthresh
            RUC37_struc(n)%ruc37(t)%acsnow        = tmp_acsnow       
            RUC37_struc(n)%ruc37(t)%sfcevp        = tmp_sfcevp       
    
            ! save output variables from local variables to global variables
            RUC37_struc(n)%ruc37(t)%shdfac        = tmp_shdfac
            RUC37_struc(n)%ruc37(t)%shdmin        = tmp_shdmin
            RUC37_struc(n)%ruc37(t)%shdmax        = tmp_shdmax
            RUC37_struc(n)%ruc37(t)%qmax          = tmp_qmax
            RUC37_struc(n)%ruc37(t)%qmin          = tmp_qmin
            RUC37_struc(n)%ruc37(t)%psis          = tmp_psis
            RUC37_struc(n)%ruc37(t)%ksat          = tmp_ksat
            RUC37_struc(n)%ruc37(t)%bclh          = tmp_bclh
            RUC37_struc(n)%ruc37(t)%qwrtz         = tmp_qwrtz
            RUC37_struc(n)%ruc37(t)%wilt          = tmp_wilt
            RUC37_struc(n)%ruc37(t)%ref           = tmp_ref
            RUC37_struc(n)%ruc37(t)%qh            = tmp_qh           
            RUC37_struc(n)%ruc37(t)%qle           = tmp_qle          
            RUC37_struc(n)%ruc37(t)%eta           = tmp_eta          
            RUC37_struc(n)%ruc37(t)%qg            = tmp_qg           
            RUC37_struc(n)%ruc37(t)%rhosnf        = tmp_rhosnf       
            RUC37_struc(n)%ruc37(t)%precipfr      = tmp_precipfr     
            RUC37_struc(n)%ruc37(t)%snomlt        = tmp_snomlt       
            RUC37_struc(n)%ruc37(t)%dew           = tmp_dew          
            RUC37_struc(n)%ruc37(t)%drip          = tmp_drip         
            RUC37_struc(n)%ruc37(t)%qs            = tmp_qs           
            RUC37_struc(n)%ruc37(t)%qsb           = tmp_qsb          
            RUC37_struc(n)%ruc37(t)%snflx         = tmp_snflx        
            RUC37_struc(n)%ruc37(t)%edir          = tmp_edir         
            RUC37_struc(n)%ruc37(t)%ec            = tmp_ec           
            RUC37_struc(n)%ruc37(t)%ett           = tmp_ett          
            RUC37_struc(n)%ruc37(t)%esnow         = tmp_esnow        
 
            
            ![ 1] output variable: emiss (unit=-). ***  surface emissivity (0.0 - 1.0).
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EMISSFORC, value = RUC37_struc(n)%ruc37(t)%emiss, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 2] output variable: ch (unit=m s-1). ***  exchange coefficient for head and moisture (m s-1).
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CH, value = RUC37_struc(n)%ruc37(t)%ch, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 3] output variable: cm (unit=m s-1). ***  exchange coefficient for momentum (m s-1).
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CM, value = RUC37_struc(n)%ruc37(t)%cm, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 4] output variable: sneqv (unit=m). ***  water equivalent of accumulated snow depth (m).
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWE, value = RUC37_struc(n)%ruc37(t)%sneqv, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            

            ![ 5] output variable: snowh (unit=m). ***  physical snow depth (m).
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDEPTH, value = RUC37_struc(n)%ruc37(t)%snowh, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 6] output variable: snowc (unit=-). ***  fractional snow cover ( fraction [0.0-1.0] )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWCOVER, value = RUC37_struc(n)%ruc37(t)%snowc, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 7] output variable: canwat (unit=kg m-2). ***  canopy moisture content (kg m-2)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPINT, value = RUC37_struc(n)%ruc37(t)%canwat, &
                                              vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 8] output variable: alb (unit=-). ***  surface albedo including possible snow-cover effect.  this is set in lsmruc,
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO, value = RUC37_struc(n)%ruc37(t)%alb, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 9] output variable: smc (unit=m^3 m-3). ***  total soil moisture content (m3 m-3)
            do i=1, RUC37_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = RUC37_struc(n)%ruc37(t)%smc(i), &
                                                  vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 10] output variable: sho (unit=m^3 m-3). ***  liquid soil moisture content (m3 m-3)
            do i=1, RUC37_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMLIQFRAC, value = RUC37_struc(n)%ruc37(t)%sho(i), &
                                                  vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 11] output variable: stc (unit=K). ***  soil temperature (k)
            do i=1, RUC37_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILTEMP, value = RUC37_struc(n)%ruc37(t)%stc(i), &
                                                  vlevel=i, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 12] output variable: tskin (unit=K). ***  skin temperature (k)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AVGSURFT, value = RUC37_struc(n)%ruc37(t)%tskin, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 13] output variable: qvg (unit=kg kg-1). ***  effective mixing ratio at the surface ( kg kg{-1} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QVG, value = RUC37_struc(n)%ruc37(t)%qvg, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 14] output variable: qcg (unit=kg kg-1). ***  effective cloud water mixing ratio at the surface ( kg kg{-1} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QCG, value = RUC37_struc(n)%ruc37(t)%qcg, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 15] output variable: qsg (unit=kg kg-1). ***  surface water vapor mixing ratio at saturation (kg kg-1)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSG, value = RUC37_struc(n)%ruc37(t)%qsg, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 16] output variable: snt75cm (unit=K). ***  snow temperature at 7.5 cm depth (k)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNNOT75CM, value = RUC37_struc(n)%ruc37(t)%snt75cm, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 17] output variable: tsnav (unit=K). ***  average snow temperature in k
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWT, value = RUC37_struc(n)%ruc37(t)%tsnav, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 18] output variable: soilm (unit=m). ***  total soil column moisture content, frozen and unfrozen ( m )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILWET, value = RUC37_struc(n)%ruc37(t)%soilm, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 19] output variable: smroot (unit=m^3 m-3). ***  available soil moisture in the root zone ( fraction [smcwlt-smcmax]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ROOTMOIST, value = RUC37_struc(n)%ruc37(t)%smroot, &
                                              vlevel=1, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 20] output variable: qh (unit=W m-2). ***  sensible heat flux ( w m{-2} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QH, value = RUC37_struc(n)%ruc37(t)%qh, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 21] output variable: qle (unit=W m-2). ***  latent heat flux (evapotranspiration) ( w m{-2} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QLE, value = RUC37_struc(n)%ruc37(t)%qle, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 22] output variable: eta (unit=kg m-2 s-1). ***  latent heat flux (evapotranspiration) ( kg m{-2} s{-1} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVAP, value = RUC37_struc(n)%ruc37(t)%eta, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 23] output variable: qg (unit=W m-2). ***  soil heat flux ( w m{-2} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QG, value = - RUC37_struc(n)%ruc37(t)%qg, &
                                              vlevel=1, unit="W m-2", direction="DN", surface_type = LIS_rc%lsm_index)
            
            ![ 24] output variable: rhosnf (unit=kg m-3). ***  density of frozen precipitation (kg m{-3})
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_FRZPREC_DEN, value = RUC37_struc(n)%ruc37(t)%rhosnf, &
                                              vlevel=1, unit="kg m-3", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 25] output variable: precipfr (unit=kg m-2). ***  time-step frozen precipitation (kg m{-2})
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_FRZPREC, value = RUC37_struc(n)%ruc37(t)%precipfr, &
                                              vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 26] output variable: snowfallac (unit=kg m-2). ***  run total snowfall accumulation (kg m{-2})
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ACC_SNOWF, value = RUC37_struc(n)%ruc37(t)%snowfallac, &
                                              vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 27] output variable: acsnow (unit=kg m-2). ***  run total frozen precipitation (kg m{-2})
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ACC_FRZPREC, value = RUC37_struc(n)%ruc37(t)%acsnow, &
                                              vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 28] output variable: sfcevp (unit=kg m-2). ***  run total evaporation flux  (kg m{-2})
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ACC_EVAP, value = RUC37_struc(n)%ruc37(t)%sfcevp, &
                                              vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 29] output variable: snomlt (unit=m). ***  snow melt water ( m )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSM, value = RUC37_struc(n)%ruc37(t)%snomlt, &
                                              vlevel=1, unit="m", direction="S2L", surface_type = LIS_rc%lsm_index)
            
            ![ 30] output variable: qsfc (unit=kg kg-1). ***  specific humidity at the surface [kg kg-1]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSFC, value = RUC37_struc(n)%ruc37(t)%qsfc, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 31] output variable: dew (unit=m). ***  dewfall (or frostfall for t<273.15) ( m )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_DEW_FROST, value = RUC37_struc(n)%ruc37(t)%dew, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 32] output variable: drip (unit=kg m-2 s-1). ***  throughfall of precipitation from canopy (kg m{-2} s{-1})
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_DRIP, value = RUC37_struc(n)%ruc37(t)%drip, &
                                              vlevel=1, unit="kg m-2 s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 33] output variable: qs (unit=kg m-2 s-1). ***  surface runoff, not infiltrating the soil ( m s{-1} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QS, value = RUC37_struc(n)%ruc37(t)%qs, &
                                              vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
            
            ![ 34] output variable: qsb (unit=kg m-2 s-1). ***  subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSB, value = RUC37_struc(n)%ruc37(t)%qsb, &
                                              vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
            
            ![ 35] output variable: snflx (unit=W m-2). ***  snow heat flux (W m-2: negative, if downward from surface)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_QH_SNOW, value = RUC37_struc(n)%ruc37(t)%snflx, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 36] output variable: edir (unit=W m-2). ***  latent heat flux component: direct soil evaporation ( w m{-2} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ESOIL, value = RUC37_struc(n)%ruc37(t)%edir, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 37] output variable: ec (unit=W m-2). ***  latent heat flux component: canopy water evaporation ( w m{-2} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ECANOP, value = RUC37_struc(n)%ruc37(t)%ec, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 38] output variable: ett (unit=W m-2). ***  latent heat flux component: total plant transpiration ( w m{-2} )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_TVEG, value = RUC37_struc(n)%ruc37(t)%ett, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 39] output variable: esnow (unit=W m-2). ***  sublimation from snowpack (w m{-2})
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SUBSNOW, value = RUC37_struc(n)%ruc37(t)%esnow, &
                                              vlevel=1, unit="W m-2", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 40] output variable: snthresh (unit=m). ***  snow depth threshold ( m )
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWTHRESH, value = RUC37_struc(n)%ruc37(t)%snthresh, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index) 
            ![ 41] output variable: shdfac (fraction). ***  shading factor for the cirrent day (fraction)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GREENNESS, value = RUC37_struc(n)%ruc37(t)%shdfac, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)


            ! reset forcing variables to zeros
            RUC37_struc(n)%ruc37(t)%lwdown = 0.0
            RUC37_struc(n)%ruc37(t)%swdown = 0.0
            RUC37_struc(n)%ruc37(t)%psurf = 0.0
            RUC37_struc(n)%ruc37(t)%rainf = 0.0
            RUC37_struc(n)%ruc37(t)%snowf = 0.0
            RUC37_struc(n)%ruc37(t)%tair = 0.0
            RUC37_struc(n)%ruc37(t)%qair = 0.0
            RUC37_struc(n)%ruc37(t)%wind_e = 0.0
            RUC37_struc(n)%ruc37(t)%wind_n = 0.0
        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        RUC37_struc(n)%forc_count = 0 
    endif ! end of alarmCheck loop 
    
    deallocate( tmp_soil_layer_thickness )
    deallocate( tmp_albedo_monthly )
    deallocate( tmp_shdfac_monthly )
    deallocate( tmp_z0brd_monthly )
    deallocate( tmp_lai_monthly )
    deallocate( tmp_smc )
    deallocate( tmp_sho )
    deallocate( tmp_stc )
    deallocate( tmp_smfr )
    deallocate( tmp_keepfr )
end subroutine RUC37_main
