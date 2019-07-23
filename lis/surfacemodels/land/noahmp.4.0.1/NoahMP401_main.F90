!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: NoahMP401_main
! \label{NoahMP401_main}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit
!  developed by Shugong Wang for the NASA Land Information System V7.
!  The initial specification of the subroutine is by Sujay Kumar.
!
!   10/25/18: Shugong Wang, Zhuo Wang; initial implementation for NoahMP401 with LIS-7
!   05/15/19: Yeosang Yoon; code added for snow DA to work
!
! !INTERFACE:
subroutine NoahMP401_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_constantsMod,  only : LIS_CONST_RHOFW   !New
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod
    use NoahMP401_lsmMod

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
!  This is the entry point for calling the NoahMP401 physics.
!  This routine calls the {\tt noahmp_driver_401} routine that performs
!  the land surface computations, to solve water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for NoahMP401
    real                 :: tmp_latitude           ! latitude in decimal degree [rad]
    real                 :: tmp_longitude          ! longitude in decimal year [rad]
    integer              :: tmp_year               ! year of the current time step [-]
    integer              :: tmp_month              ! month of the current time step [-]
    integer              :: tmp_day                ! day of the current time step [-]
    integer              :: tmp_hour               ! hour of the current time step [-]
    integer              :: tmp_minute             ! minute of the current time step [-]
    real                 :: tmp_dz8w               ! reference height of temperature and humidity [m]
    real                 :: tmp_dt                 ! timestep [s]
    integer              :: tmp_ttile              ! tile No. [-]
    integer              :: tmp_itimestep          ! timestep number [-]

    real, allocatable    :: tmp_sldpth(:)          ! thickness of soil layers [m]
    integer              :: tmp_nsoil              ! number of soil layers [-]
    integer              :: tmp_nsnow              ! maximum number of snow layers (e.g. 3) [-]
    integer              :: tmp_vegetype           ! vegetation type [-]
    integer              :: tmp_soiltype           ! soil type [-]
    real, allocatable    :: tmp_shdfac_monthly(:)  ! monthly values for green vegetation fraction []
    real                 :: tmp_tbot               ! deep soil temperature [K]
    integer              :: tmp_urban_vegetype     ! urban land cover type index [-]
    integer              :: tmp_cropcat            ! crop category [-]
    real                 :: tmp_planting           ! planting date [-]
    real                 :: tmp_harvest            ! harvest date [-]
    real                 :: tmp_season_gdd         ! growing season GDD [-]
    character            :: tmp_landuse_tbl_name   ! Noah model landuse parameter table [-]
    character            :: tmp_soil_tbl_name      ! Noah model soil parameter table [-]
    character            :: tmp_gen_tbl_name       ! Noah model general parameter table [-]
    character            :: tmp_noahmp_tbl_name    ! NoahMP parameter table [-]
    integer              :: tmp_dveg_opt           ! dynamic vegetation, (1->off; 2->on); with opt_crs=1 [-]
    integer              :: tmp_crs_opt            ! canopt stomatal resistance (1->Ball-Berry; 2->Jarvis) [-]
    integer              :: tmp_btr_opt            ! soil moisture factor for stomatal resistance (1->Noah;2->CLM;3->SSiB) [-]
    integer              :: tmp_run_opt            ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS) [-]
    integer              :: tmp_sfc_opt            ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97) [-]
    integer              :: tmp_frz_opt            ! supercooled liquid water (1->NY06; 2->Koren99) [-]
    integer              :: tmp_inf_opt            ! frozen soil permeability (1->NY06; 2->Koren99) [-]
    integer              :: tmp_rad_opt            ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg) [-]
    integer              :: tmp_alb_opt            ! snow surface albedo (1->BATS; 2->CLASS) [-]
    integer              :: tmp_snf_opt            ! rainfall & snowfall (1->Jordan91; 2->BATS; 3->Noah) [-]
    integer              :: tmp_tbot_opt           ! lower boundary of soil temperature [-]
    integer              :: tmp_stc_opt            ! snow/soil temperature time scheme [-]
    integer              :: tmp_gla_opt            ! glacier option (1->phase change; 2->simple) [-]
    integer              :: tmp_rsf_opt            ! surface resistance (1->Sakaguchi/Zeng;2->Seller;3->mod Sellers;4->1+snow) [-]
    integer              :: tmp_soil_opt           ! soil configuration option [-]
    integer              :: tmp_pedo_opt           ! soil pedotransfer function option [-]
    integer              :: tmp_crop_opt           ! crop model option (0->none; 1->Liu et al.; 2->Gecros) [-]
    integer              :: tmp_iz0tlnd            ! option of Chen adjustment of Czil (not used) [-]
    integer              :: tmp_urban_opt          ! urban physics option [-]
    real, allocatable    :: tmp_soilcomp(:)        ! soil sand and clay percentage [-]
    real                 :: tmp_soilcL1            ! soil texture in layer 1 [-]
    real                 :: tmp_soilcL2            ! soil texture in layer 2 [-]
    real                 :: tmp_soilcL3            ! soil texture in layer 3 [-]
    real                 :: tmp_soilcL4            ! soil texture in layer 4 [-]
    real                 :: tmp_tair               ! air temperature [K]
    real                 :: tmp_psurf              ! air pressure [Pa]
    real                 :: tmp_wind_e             ! U wind component [m/s]
    real                 :: tmp_wind_n             ! V wind component [m/s]
    real                 :: tmp_qair               ! specific humidity [kg/kg]
    real                 :: tmp_swdown             ! downward solar radiation [W m-2]
    real                 :: tmp_lwdown             ! downward longwave radiation [W m-2]
    real                 :: tmp_prcp               ! total precipitation (rainfall+snowfall) [mm]
    real                 :: tmp_tsk                ! surface radiative temperature [K]
    real                 :: tmp_hfx                ! sensible heat flux [W m-2]
!   real                 :: tmp_fsh                ! sensible heat flux [W/m2]

    real                 :: tmp_qfx                ! latent heat flux [kg s-1 m-2]
    real                 :: tmp_lh                 ! latent heat flux [W m-2]
    real                 :: tmp_grdflx             ! ground/snow heat flux [W m-2]
    real                 :: tmp_sfcrunoff          ! accumulated surface runoff [m]
    real                 :: tmp_udrrunoff          ! accumulated sub-surface runoff [m]
    real                 :: tmp_albedo             ! total grid albedo [-]

    real                 :: tmp_snowc              ! snow cover fraction [-]
    real, allocatable    :: tmp_smc(:)             ! volumetric soil moisture [m3/m3]
    real, allocatable    :: tmp_sh2o(:)            ! volumetric liquid soil moisture [m3/m3]
    real, allocatable    :: tmp_tslb(:)            ! soil temperature [K]
    real                 :: tmp_sneqv              ! snow water equivalent [mm]
    real                 :: tmp_snowh              ! physical snow depth [m]
    real                 :: tmp_canwat             ! total canopy water + ice [mm]
    real                 :: tmp_acsnom             ! accumulated snow melt leaving pack [-]
    real                 :: tmp_acsnow             ! accumulated snow on grid [mm]
    real                 :: tmp_emiss              ! surface bulk emissivity [-]
    real                 :: tmp_rs                 ! total stomatal resistance [s/m]
    integer              :: tmp_isnow              ! actual no. of snow layers [-]
    real                 :: tmp_tv                 ! vegetation leaf temperature [K]
    real                 :: tmp_tg                 ! bulk ground surface temperature [K]
    real                 :: tmp_canice             ! canopy-intercepted ice [mm]
    real                 :: tmp_canliq             ! canopy-intercepted liquid water [mm]
    real                 :: tmp_eah                ! canopy air vapor pressure [Pa]
    real                 :: tmp_tah                ! canopy air temperature [K]
    real                 :: tmp_cm                 ! bulk momentum drag coefficient [-]
    real                 :: tmp_ch                 ! bulk sensible heat exchange coefficient [-]
    real                 :: tmp_fwet               ! wetted or snowed fraction of canopy [-]
    real                 :: tmp_sneqvo             ! snow mass at last time step [mm h2o]
    real                 :: tmp_albold             ! snow albedo at last time step [-]
    real                 :: tmp_qsnow              ! snowfall on the ground [mm/s]
    real                 :: tmp_wslake             ! lake water storage [mm]
    real                 :: tmp_zwt                ! water table depth [m]
    real                 :: tmp_wa                 ! water in the "aquifer" [mm]
    real                 :: tmp_wt                 ! water in aquifer and saturated soil [mm]
    real, allocatable    :: tmp_tsno(:)            ! snow layer temperature [K]
    real, allocatable    :: tmp_zss(:)             ! snow/soil layer depth from snow surface [m]
    real, allocatable    :: tmp_snowice(:)         ! snow layer ice [mm]
    real, allocatable    :: tmp_snowliq(:)         ! snow layer liquid water [mm]
    real                 :: tmp_lfmass             ! leaf mass [g/m2]
    real                 :: tmp_rtmass             ! mass of fine roots [g/m2]
    real                 :: tmp_stmass             ! stem mass [g/m2]
    real                 :: tmp_wood               ! mass of wood (including woody roots) [g/m2]
    real                 :: tmp_stblcp             ! stable carbon in deep soil [g/m2]
    real                 :: tmp_fastcp             ! short-lived carbon in shallow soil [g/m2]
    real                 :: tmp_lai                ! leaf area index [-]
    real                 :: tmp_sai                ! stem area index [-]
    real                 :: tmp_tauss              ! snow age factor [-]
    real, allocatable    :: tmp_smoiseq(:)         ! equilibrium volumetric soil moisture content [m3/m3]
    real                 :: tmp_smcwtd             ! soil moisture content in the layer to the water table when deep [-]
    real                 :: tmp_deeprech           ! recharge to the water table when deep [-]
    real                 :: tmp_rech               ! recharge to the water table (diagnostic) [-]
    real                 :: tmp_grain              ! mass of grain XING [g/m2]
    real                 :: tmp_gdd                ! growing degree days XING (based on 10C) [-]
    integer              :: tmp_pgs                ! growing degree days XING [-]
    real, allocatable    :: tmp_gecros_state(:)    ! optional gecros crop [-]
    real                 :: tmp_t2mv               ! 2m temperature of vegetation part [K]
    real                 :: tmp_t2mb               ! 2m temperature of bare ground part [K]
    real                 :: tmp_q2mv               ! 2m mixing ratio of vegetation part [-]
    real                 :: tmp_q2mb               ! 2m mixing ratio of bare ground part [-]
    real                 :: tmp_trad               ! surface radiative temperature [K]
    real                 :: tmp_nee                ! net ecosys exchange of CO2 [g/m2/s CO2]
    real                 :: tmp_gpp                ! gross primary assimilation of carbon [g/m2/s C]
    real                 :: tmp_npp                ! net primary productivity of carbon [g/m2/s C]
    real                 :: tmp_fveg               ! Noah-MP green vegetation fraction [-]
    real                 :: tmp_runsf              ! surface runoff [mm/s]
    real                 :: tmp_runsb              ! subsurface runoff [mm/s]
    real                 :: tmp_ecan               ! evaporation of intercepted water [mm/s]
    real                 :: tmp_edir               ! soil surface evaporation rate [mm/s]
    real                 :: tmp_etran              ! transpiration rate [mm/s]
    real                 :: tmp_rainf              ! rainfall rate [kg s-1]
    real                 :: tmp_snowf              ! snowfall rate [kg s-1]
    real                 :: tmp_fsa                ! total absorbed solar radiation [W/m2]
    real                 :: tmp_fira               ! total net longwave radiation [+ to atm] [W/m2]
    real                 :: tmp_apar               ! photosyn active energy by canopy [W/m2]
    real                 :: tmp_psn                ! total photosynthesis [+] [umol co2/m2/s]
    real                 :: tmp_sav                ! solar radiation absorbed by vegetation [W/m2]
    real                 :: tmp_sag                ! solar radiation absorbed by ground [W/m2]
    real                 :: tmp_rssun              ! sunlit leaf stomatal resistance [s/m]
    real                 :: tmp_rssha              ! shaded leaf stomatal resistance [s/m]
    real                 :: tmp_bgap               ! between gap fraction [-]
    real                 :: tmp_wgap               ! within gap fraction [-]
    real                 :: tmp_tgb                ! bare ground temperature [K]
    real                 :: tmp_tgv                ! under canopy ground temperature [K]
    real                 :: tmp_chv                ! sensible heat exchange coefficient vegetated [-]
    real                 :: tmp_chb                ! sensible heat exchange coefficient bare-ground [-]
    real                 :: tmp_shg                ! veg ground sensible heat [+ to atm] [W/m2]
    real                 :: tmp_shc                ! canopy sensible heat [+ to atm] [W/m2]
    real                 :: tmp_shb                ! bare sensible heat [+ to atm] [W/m2]
    real                 :: tmp_evg                ! veg ground evaporation [+ to atm] [W/m2]
    real                 :: tmp_evb                ! bare soil evaporation [+ to atm] [W/m2]
    real                 :: tmp_ghv                ! veg ground heat flux [+ to soil] [W/m2]
    real                 :: tmp_ghb                ! bare ground heat flux [+ to soil] [W/m2]
    real                 :: tmp_irg                ! veg ground net LW radiation [+ to atm] [W/m2]
    real                 :: tmp_irc                ! canopy net LW radiation [+ to atm] [W/m2]
    real                 :: tmp_irb                ! bare net LW radiation [+ to atm] [W/m2]
    real                 :: tmp_tr                 ! transpiration [ to atm] [W/m2]
    real                 :: tmp_evc                ! canopy evaporation heat [to atm] [W/m2]
    real                 :: tmp_chleaf             ! leaf exchange coefficient [-]
    real                 :: tmp_chuc               ! under canopy exchange coefficient [-]
    real                 :: tmp_chv2               ! veg 2m exchange coefficient [-]
    real                 :: tmp_chb2               ! bare 2m exchange coefficient [-]
    real                 :: tmp_qsnbot             ! melting water out of snow bottom [kg m-2 s-1]
    real                 :: tmp_subsnow            ! snow sublimation rate [kg m-2 s-1]
    real                 :: tmp_pah                ! precipitation advected heat - total (W/m2)

    ! Code added by Zhuo Wang 02/28/2019
    real                 :: AvgSurfT_out           ! average surface temperature [K]
    real                 :: TWS_out                ! terrestrial water storage [mm]
    ! Code added by David Mocko 04/25/2019
    real                 :: startsm, startswe, startint, startgw, endsm

    allocate( tmp_sldpth( NOAHMP401_struc(n)%nsoil ) )
    allocate( tmp_shdfac_monthly( 12 ) )
    allocate( tmp_soilcomp( 8 ) )
    allocate( tmp_smc( NOAHMP401_struc(n)%nsoil ) )
    allocate( tmp_sh2o( NOAHMP401_struc(n)%nsoil ) )
    allocate( tmp_tslb( NOAHMP401_struc(n)%nsoil ) )
    allocate( tmp_tsno( NOAHMP401_struc(n)%nsnow ) )
    allocate( tmp_zss( NOAHMP401_struc(n)%nsnow+NOAHMP401_struc(n)%nsoil) )
    allocate( tmp_snowice( NOAHMP401_struc(n)%nsnow ) )
    allocate( tmp_snowliq( NOAHMP401_struc(n)%nsnow ) )
    allocate( tmp_smoiseq( NOAHMP401_struc(n)%nsoil ) )
    allocate( tmp_gecros_state( 60 ) )

    ! check NoahMP401 alarm. If alarm is ring, run model.

    alarmCheck = LIS_isAlarmRinging(LIS_rc, "NoahMP401 model alarm")

    if (alarmCheck) Then
        do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            dt = LIS_rc%ts
            row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
            lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

            ! retrieve forcing data from NOAHMP401_struc(n)%noahmp401(t) and assign to local variables
            ! tair: air temperature
            tmp_tair       = NOAHMP401_struc(n)%noahmp401(t)%tair   / NOAHMP401_struc(n)%forc_count
            ! Yeosang Yoon, for snow DA
            NOAHMP401_struc(n)%noahmp401(t)%sfctmp = tmp_tair

            ! psurf: air pressure
            tmp_psurf      = NOAHMP401_struc(n)%noahmp401(t)%psurf  / NOAHMP401_struc(n)%forc_count

            ! wind_e: U wind component
            tmp_wind_e     = NOAHMP401_struc(n)%noahmp401(t)%wind_e / NOAHMP401_struc(n)%forc_count

            ! wind_n: V wind component
            tmp_wind_n     = NOAHMP401_struc(n)%noahmp401(t)%wind_n / NOAHMP401_struc(n)%forc_count

            ! qair: specific humidity
            tmp_qair       = NOAHMP401_struc(n)%noahmp401(t)%qair   / NOAHMP401_struc(n)%forc_count

            ! swdown: downward solar radiation
            tmp_swdown     = NOAHMP401_struc(n)%noahmp401(t)%swdown / NOAHMP401_struc(n)%forc_count

            ! lwdown: downward longwave radiation
            tmp_lwdown     = NOAHMP401_struc(n)%noahmp401(t)%lwdown / NOAHMP401_struc(n)%forc_count

            ! prcp: total precipitation (rainfall+snowfall)
            ! Both NoahMP-3.6.1 and NoahMP-4.0.1 requires total precipitation as forcing input.
            ! In LIS/NoahMP-3.6.1, the input forcing is total precipitation [mm], but in
            ! LIS/NoahMP-4.0.1, the forcing data provides precipitation rate [mm/s] !!!
            tmp_prcp       = dt * (NOAHMP401_struc(n)%noahmp401(t)%prcp   / NOAHMP401_struc(n)%forc_count)

            ! check validity of tair
            if(tmp_tair .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable tair in NoahMP401"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of psurf
            if(tmp_psurf .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable psurf in NoahMP401"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of wind_e
            if(tmp_wind_e .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable wind_e in NoahMP401"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of wind_n
            if(tmp_wind_n .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable wind_n in NoahMP401"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of qair
            if(tmp_qair .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable qair in NoahMP401"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of swdown
            if(tmp_swdown .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable swdown in NoahMP401"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of lwdown
            if(tmp_lwdown .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable lwdown in NoahMP401"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of prcp
            if(tmp_prcp .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable prcp in NoahMP401"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            !

            tmp_latitude  = lat
            tmp_longitude = lon
            tmp_year   = LIS_rc%yr
            tmp_month  = LIS_rc%mo
            tmp_day    = LIS_rc%da
            tmp_hour   = LIS_rc%hr
            tmp_minute = LIS_rc%mn

            ! Added by Zhuo Wang on 11/13/2018
            tmp_ttile  = t
            tmp_itimestep  = LIS_rc%tscount(n)

            ! get parameters
            tmp_dt                = NOAHMP401_struc(n)%ts
            tmp_sldpth(:)         = NOAHMP401_struc(n)%sldpth(:)
            tmp_nsoil             = NOAHMP401_struc(n)%nsoil
            tmp_nsnow             = NOAHMP401_struc(n)%nsnow
            tmp_vegetype          = NOAHMP401_struc(n)%noahmp401(t)%vegetype
            tmp_soiltype          = NOAHMP401_struc(n)%noahmp401(t)%soiltype
            ! Multiply shdfac by 100.0 because noahmpdrv.f90
            ! expects it in units of percentage, not fraction.
            tmp_shdfac_monthly(:) = NOAHMP401_struc(n)%noahmp401(t)%shdfac_monthly(:) * 100.0
            tmp_tbot              = NOAHMP401_struc(n)%noahmp401(t)%tbot
            tmp_urban_vegetype    = LIS_rc%urbanclass
            tmp_cropcat           = LIS_rc%cropclass
            tmp_dveg_opt          = NOAHMP401_struc(n)%dveg_opt
            tmp_crs_opt           = NOAHMP401_struc(n)%crs_opt
            tmp_btr_opt           = NOAHMP401_struc(n)%btr_opt
            tmp_run_opt           = NOAHMP401_struc(n)%run_opt
            tmp_sfc_opt           = NOAHMP401_struc(n)%sfc_opt
            tmp_frz_opt           = NOAHMP401_struc(n)%frz_opt
            tmp_inf_opt           = NOAHMP401_struc(n)%inf_opt
            tmp_rad_opt           = NOAHMP401_struc(n)%rad_opt
            tmp_alb_opt           = NOAHMP401_struc(n)%alb_opt
            tmp_snf_opt           = NOAHMP401_struc(n)%snf_opt
            tmp_tbot_opt          = NOAHMP401_struc(n)%tbot_opt
            tmp_stc_opt           = NOAHMP401_struc(n)%stc_opt
            tmp_gla_opt           = NOAHMP401_struc(n)%gla_opt
            tmp_rsf_opt           = NOAHMP401_struc(n)%rsf_opt
            tmp_soil_opt          = NOAHMP401_struc(n)%soil_opt
            tmp_pedo_opt          = NOAHMP401_struc(n)%pedo_opt
            tmp_crop_opt          = NOAHMP401_struc(n)%crop_opt
            tmp_iz0tlnd           = 0
            tmp_urban_opt         = NOAHMP401_struc(n)%urban_opt
! Multiply reference height by 2.0 because module_sf_noahmpdrv
! expects this variable to be in terms of a thickness of the
! atmospheric layers, and it later divides this value by 2.0.
! Thus, the LIS user should specify the exact height of the
! reference in lis.config, and module_sf_noahmpdrv will then
! correctly use this actual value.  This code is confirmed in
! the HRLDAS driver, which also multiplies this value by 2.0.
! 11/30/2018 - dmm
            tmp_dz8w              = NOAHMP401_struc(n)%dz8w * 2.0

            if (tmp_crop_opt.ne.0) then 
               tmp_planting   = NOAHMP401_struc(n)%noahmp401(t)%planting
               tmp_harvest    = NOAHMP401_struc(n)%noahmp401(t)%harvest
               tmp_season_gdd = NOAHMP401_struc(n)%noahmp401(t)%season_gdd
               if (tmp_crop_opt.eq.2) then
                  tmp_gecros_state(:) = NOAHMP401_struc(n)%noahmp401(t)%gecros_state(:)
               else
                  tmp_gecros_state(:) = 0.0
               endif
            else
               tmp_planting   = 0.0
               tmp_harvest    = 0.0
               tmp_season_gdd = 0.0
            endif

! Zhuo Wang tested on 11/15/2018, not read from LDT-generated netcdf input file
            if (tmp_soil_opt.eq.2) then 
               tmp_soilcL1 = NOAHMP401_struc(n)%noahmp401(t)%soilcL1
               tmp_soilcL2 = NOAHMP401_struc(n)%noahmp401(t)%soilcL2
               tmp_soilcL3 = NOAHMP401_struc(n)%noahmp401(t)%soilcL3
               tmp_soilcL4 = NOAHMP401_struc(n)%noahmp401(t)%soilcL4
            else
               tmp_soilcL1 = 0.0
               tmp_soilcL2 = 0.0
               tmp_soilcL3 = 0.0
               tmp_soilcL4 = 0.0
            endif
            if (tmp_soil_opt.eq.3) then 
               tmp_soilcomp(:) = NOAHMP401_struc(n)%noahmp401(t)%soilcomp(:)
            else
               tmp_soilcomp = 0.0
            endif

            ! get state variables
            tmp_sfcrunoff       = NOAHMP401_struc(n)%noahmp401(t)%sfcrunoff
            tmp_udrrunoff       = NOAHMP401_struc(n)%noahmp401(t)%udrrunoff
            tmp_smc(:)          = NOAHMP401_struc(n)%noahmp401(t)%smc(:)
            tmp_sh2o(:)         = NOAHMP401_struc(n)%noahmp401(t)%sh2o(:)
            tmp_tslb(:)         = NOAHMP401_struc(n)%noahmp401(t)%tslb(:)
            tmp_sneqv           = NOAHMP401_struc(n)%noahmp401(t)%sneqv
            tmp_snowh           = NOAHMP401_struc(n)%noahmp401(t)%snowh
            tmp_canwat          = NOAHMP401_struc(n)%noahmp401(t)%canwat
            tmp_acsnom          = NOAHMP401_struc(n)%noahmp401(t)%acsnom
            tmp_acsnow          = NOAHMP401_struc(n)%noahmp401(t)%acsnow
            tmp_isnow           = NOAHMP401_struc(n)%noahmp401(t)%isnow
            tmp_tv              = NOAHMP401_struc(n)%noahmp401(t)%tv
            tmp_tg              = NOAHMP401_struc(n)%noahmp401(t)%tg
            tmp_canice          = NOAHMP401_struc(n)%noahmp401(t)%canice
            tmp_canliq          = NOAHMP401_struc(n)%noahmp401(t)%canliq
            tmp_eah             = NOAHMP401_struc(n)%noahmp401(t)%eah
            tmp_tah             = NOAHMP401_struc(n)%noahmp401(t)%tah
            tmp_cm              = NOAHMP401_struc(n)%noahmp401(t)%cm
            tmp_ch              = NOAHMP401_struc(n)%noahmp401(t)%ch
            tmp_fwet            = NOAHMP401_struc(n)%noahmp401(t)%fwet
            tmp_sneqvo          = NOAHMP401_struc(n)%noahmp401(t)%sneqvo
            tmp_albold          = NOAHMP401_struc(n)%noahmp401(t)%albold
            tmp_qsnow           = NOAHMP401_struc(n)%noahmp401(t)%qsnow
            tmp_wslake          = NOAHMP401_struc(n)%noahmp401(t)%wslake
            tmp_zwt             = NOAHMP401_struc(n)%noahmp401(t)%zwt
            tmp_wa              = NOAHMP401_struc(n)%noahmp401(t)%wa
            tmp_wt              = NOAHMP401_struc(n)%noahmp401(t)%wt
            tmp_tsno(:)         = NOAHMP401_struc(n)%noahmp401(t)%tsno(:)
            tmp_zss(:)          = NOAHMP401_struc(n)%noahmp401(t)%zss(:)
            tmp_snowice(:)      = NOAHMP401_struc(n)%noahmp401(t)%snowice(:)
            tmp_snowliq(:)      = NOAHMP401_struc(n)%noahmp401(t)%snowliq(:)
            tmp_lfmass          = NOAHMP401_struc(n)%noahmp401(t)%lfmass
            tmp_rtmass          = NOAHMP401_struc(n)%noahmp401(t)%rtmass
            tmp_stmass          = NOAHMP401_struc(n)%noahmp401(t)%stmass
            tmp_wood            = NOAHMP401_struc(n)%noahmp401(t)%wood
            tmp_stblcp          = NOAHMP401_struc(n)%noahmp401(t)%stblcp
            tmp_fastcp          = NOAHMP401_struc(n)%noahmp401(t)%fastcp
            tmp_lai             = NOAHMP401_struc(n)%noahmp401(t)%lai
            tmp_sai             = NOAHMP401_struc(n)%noahmp401(t)%sai
            tmp_tauss           = NOAHMP401_struc(n)%noahmp401(t)%tauss
            tmp_smoiseq(:)      = NOAHMP401_struc(n)%noahmp401(t)%smoiseq(:)
            tmp_smcwtd          = NOAHMP401_struc(n)%noahmp401(t)%smcwtd
            tmp_deeprech        = NOAHMP401_struc(n)%noahmp401(t)%deeprech
            tmp_rech            = NOAHMP401_struc(n)%noahmp401(t)%rech
            tmp_grain           = NOAHMP401_struc(n)%noahmp401(t)%grain
            tmp_gdd             = NOAHMP401_struc(n)%noahmp401(t)%gdd
            tmp_pgs             = NOAHMP401_struc(n)%noahmp401(t)%pgs

! Calculate water storages at start of timestep
            startsm = 0.0
            do i = 1,tmp_nsoil
               startsm = startsm +                                     &
                          (tmp_smc(i) * tmp_sldpth(i) * LIS_CONST_RHOFW)
            enddo
            startswe = tmp_sneqv
            startint = tmp_canliq + tmp_canice
            startgw  = tmp_wa

            ! call model physics

            call noahmp_driver_401(n                     , & ! in    - nest id [-]
                                   tmp_ttile             , & ! in    - tile id [-]
                                   tmp_itimestep         , & ! in    - timestep number [-]
                                   tmp_latitude          , & ! in    - latitude in decimal degree [rad]
                                   tmp_longitude         , & ! in    - longitude in decimal year [rad]
                                   tmp_year              , & ! in    - year of the currrent time step [-]
                                   tmp_month             , & ! in    - month of the current time step [-]
                                   tmp_day               , & ! in    - day of the current time step [-]
                                   tmp_hour              , & ! in    - hour of the current time step [-]
                                   tmp_minute            , & ! in    - minute of the current time step [-]
                                   tmp_dz8w              , & ! in    - thickness of atmospheric layers [m]
                                   tmp_dt                , & ! in    - timestep [s]
                                   tmp_sldpth            , & ! in    - thickness of soil layers [m]
                                   tmp_nsoil             , & ! in    - number of soil layers [-]
                                   tmp_nsnow             , & ! in    - maximum number of snow layers (e.g. 3) [-]
                                   tmp_vegetype          , & ! in    - vegetation type [-]
                                   tmp_soiltype          , & ! in    - soil type [-]
                                   tmp_shdfac_monthly    , & ! in    - monthly values for green vegetation fraction []
                                   tmp_tbot              , & ! in    - deep soil temperature [K]
                                   tmp_urban_vegetype    , & ! in    - urban land cover type index [-]
                                   tmp_cropcat           , & ! in    - crop category [-]
                                   tmp_planting          , & ! in    - planting date [-]
                                   tmp_harvest           , & ! in    - harvest date [-]
                                   tmp_season_gdd        , & ! in    - growing season GDD [-]
                                   tmp_dveg_opt          , & ! in    - dynamic vegetation, 1->off; 2->on); with opt_crs=1 [-]
                                   tmp_crs_opt           , & ! in    - canopt stomatal resistance (1->Ball-Berry; 2->Jarvis) [-]
                                   tmp_btr_opt           , & ! in    - soil moisture factor for stomatal resistance(1->Noah;2->CLM;3->SSiB) [-]
                                   tmp_run_opt           , & ! in    - runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS) [-]
                                   tmp_sfc_opt           , & ! in    - surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97) [-]
                                   tmp_frz_opt           , & ! in    - supercooled liquid water (1->NY06; 2->Koren99) [-]
                                   tmp_inf_opt           , & ! in    - frozen soil permeability (1->NY06; 2->Koren99) [-]
                                   tmp_rad_opt           , & ! in    - radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg) [-]
                                   tmp_alb_opt           , & ! in    - snow surface albedo (1->BATS; 2->CLASS) [-]
                                   tmp_snf_opt           , & ! in    - rainfall & snowfall (1->Jordan91; 2->BATS; 3->Noah) [-]
                                   tmp_tbot_opt          , & ! in    - lower boundary of soil temperature [-]
                                   tmp_stc_opt           , & ! in    - snow/soil temperature time scheme [-]
                                   tmp_gla_opt           , & ! in    - glacier option (1->phase change; 2->simple) [-]
                                   tmp_rsf_opt           , & ! in    - surface resistance(1->Sakaguchi/Zeng;2->Seller;3->mod Sellers;4->1+snow) [-]
                                   tmp_soil_opt          , & ! in    - soil configuration option [-]
                                   tmp_pedo_opt          , & ! in    - soil pedotransfer function option [-]
                                   tmp_crop_opt          , & ! in    - crop model option (0->none; 1->Liu et al.; 2->Gecros) [-]
                                   tmp_iz0tlnd           , & ! in    - option of Chen adjustment of Czil (not used) [-]
                                   tmp_urban_opt         , & ! in    - urban physics option [-]
                                   tmp_soilcomp          , & ! in    - soil sand and clay percentage [-]
                                   tmp_soilcL1           , & ! in    - soil texture in layer 1 [-]
                                   tmp_soilcL2           , & ! in    - soil texture in layer 2 [-]
                                   tmp_soilcL3           , & ! in    - soil texture in layer 3 [-]
                                   tmp_soilcL4           , & ! in    - soil texture in layer 4 [-]
                                   tmp_tair              , & ! in    - air temperature [K]
                                   tmp_psurf             , & ! in    - air pressure [Pa]
                                   tmp_wind_e            , & ! in    - U wind component [m/s]
                                   tmp_wind_n            , & ! in    - V wind component [m/s]
                                   tmp_qair              , & ! in    - specific humidity [kg/kg]
                                   tmp_swdown            , & ! in    - downward solar radiation [W m-2]
                                   tmp_lwdown            , & ! in    - downward longwave radiation [W m-2]
                                   tmp_prcp              , & ! in    - total precipitation (rainfall+snowfall) [mm]
                                   tmp_tsk               , & ! out   - surface radiative temperature [K]
                                   tmp_hfx               , & ! out   - sensible heat flux [W m-2]
!                                  tmp_fsh               , & ! out   - sensible heat flux [W/m2]

                                   tmp_qfx               , & ! out   - latent heat flux [kg s-1 m-2]
                                   tmp_lh                , & ! out   - latent heat flux [W m-2]
                                   tmp_grdflx            , & ! out   - ground/snow heat flux [W m-2]
                                   tmp_sfcrunoff         , & ! inout - accumulated surface runoff [m]
                                   tmp_udrrunoff         , & ! inout - accumulated sub-surface runoff [m]
                                   tmp_albedo            , & ! out   - total grid albedo [-]
                                   tmp_qsnbot            , & ! out   - melting water out of snow bottom [kg m-2 s-1]
                                   tmp_subsnow           , & ! out   - snow sublimation rate [kg m-2 s-1]
                                   tmp_snowc             , & ! out   - snow cover fraction [-]
                                   tmp_smc               , & ! inout - volumetric soil moisture [m3/m3]
                                   tmp_pah               , & ! out   - precipitation advected heat - total (W/m2)
                                   tmp_sh2o              , & ! inout - volumetric liquid soil moisture [m3/m3]
                                   tmp_tslb              , & ! inout - soil temperature [K]
                                   tmp_sneqv             , & ! inout - snow water equivalent [mm]
                                   tmp_snowh             , & ! inout - physical snow depth [m]
                                   tmp_canwat            , & ! inout - total canopy water + ice [mm]
                                   tmp_acsnom            , & ! inout - accumulated snow melt leaving pack [-]
                                   tmp_acsnow            , & ! inout - accumulated snow on grid [mm]
                                   tmp_emiss             , & ! out   - surface bulk emissivity [-]
                                   tmp_rs                , & ! out   - total stomatal resistance [s/m]
                                   tmp_isnow             , & ! inout - actual no. of snow layers [-]
                                   tmp_tv                , & ! inout - vegetation leaf temperature [K]
                                   tmp_tg                , & ! inout - bulk ground surface temperature [K]
                                   tmp_canice            , & ! inout - canopy-intercepted ice [mm]
                                   tmp_canliq            , & ! inout - canopy-intercepted liquid water [mm]
                                   tmp_eah               , & ! inout - canopy air vapor pressure [Pa]
                                   tmp_tah               , & ! inout - canopy air temperature [K]
                                   tmp_cm                , & ! inout - bulk momentum drag coefficient [-]
                                   tmp_ch                , & ! inout - bulk sensible heat exchange coefficient [-]
                                   tmp_fwet              , & ! inout - wetted or snowed fraction of canopy [-]
                                   tmp_sneqvo            , & ! inout - snow mass at last time step [mm h2o]
                                   tmp_albold            , & ! inout - snow albedo at last time step [-]
                                   tmp_qsnow             , & ! inout - snowfall on the ground [mm/s]
                                   tmp_wslake            , & ! inout - lake water storage [mm]
                                   tmp_zwt               , & ! inout - water table depth [m]
                                   tmp_wa                , & ! inout - water in the "aquifer" [mm]
                                   tmp_wt                , & ! inout - water in aquifer and saturated soil [mm]
                                   tmp_tsno              , & ! inout - snow layer temperature [K]
                                   tmp_zss               , & ! inout - snow/soil layer depth from snow surface [m]
                                   tmp_snowice           , & ! inout - snow layer ice [mm]
                                   tmp_snowliq           , & ! inout - snow layer liquid water [mm]
                                   tmp_lfmass            , & ! inout - leaf mass [g/m2]
                                   tmp_rtmass            , & ! inout - mass of fine roots [g/m2]
                                   tmp_stmass            , & ! inout - stem mass [g/m2]
                                   tmp_wood              , & ! inout - mass of wood (including woody roots) [g/m2]
                                   tmp_stblcp            , & ! inout - stable carbon in deep soil [g/m2]
                                   tmp_fastcp            , & ! inout - short-lived carbon in shallow soil [g/m2]
                                   tmp_lai               , & ! inout - leaf area index [-]
                                   tmp_sai               , & ! inout - stem area index [-]
                                   tmp_tauss             , & ! inout - snow age factor [-]
                                   tmp_smoiseq           , & ! inout - equilibrium volumetric soil moisture content [m3/m3]
                                   tmp_smcwtd            , & ! inout - soil moisture content in the layer to the water table when deep [-]
                                   tmp_deeprech          , & ! inout - recharge to the water table when deep [-]
                                   tmp_rech              , & ! inout - recharge to the water table (diagnostic) [-]
                                   tmp_grain             , & ! inout - mass of grain XING [g/m2]
                                   tmp_gdd               , & ! inout - growing degree days XING (based on 10C) [-]
                                   tmp_pgs               , & ! inout - growing degree days XING [-]
                                   tmp_gecros_state      , & ! inout - optional gecros crop [-]
                                   tmp_t2mv              , & ! out   - 2m temperature of vegetation part [K]
                                   tmp_t2mb              , & ! out   - 2m temperature of bare ground part [K]
                                   tmp_q2mv              , & ! out   - 2m mixing ratio of vegetation part [-]
                                   tmp_q2mb              , & ! out   - 2m mixing ratio of bare ground part [-]
                                   tmp_trad              , & ! out   - surface radiative temperature [K]
                                   tmp_nee               , & ! out   - net ecosys exchange of CO2 [g/m2/s CO2]
                                   tmp_gpp               , & ! out   - gross primary assimilation of carbon [g/m2/s C]
                                   tmp_npp               , & ! out   - net primary productivity of carbon [g/m2/s C]
                                   tmp_fveg              , & ! out   - Noah-MP green vegetation fraction [-]
                                   tmp_runsf             , & ! out   - surface runoff [mm/s]
                                   tmp_runsb             , & ! out   - subsurface runoff [mm/s]
                                   tmp_ecan              , & ! out   - evaporation of intercepted water [mm/s]
                                   tmp_edir              , & ! out   - soil surface evaporation rate [mm/s]
                                   tmp_etran             , & ! out   - transpiration rate [mm/s]
                                   tmp_rainf             , & ! out   - raifall rate [km s-1]
                                   tmp_snowf             , & ! out   - snowfall rate [kg s-1]
                                   tmp_fsa               , & ! out   - total absorbed solar radiation [W/m2]
                                   tmp_fira              , & ! out   - total net longwave radiation [+ to atm] [W/m2]
                                   tmp_apar              , & ! out   - photosyn active energy by canopy [W/m2]
                                   tmp_psn               , & ! out   - total photosynthesis [+] [umol co2/m2/s]
                                   tmp_sav               , & ! out   - solar radiation absorbed by vegetation [W/m2]
                                   tmp_sag               , & ! out   - solar radiatiob absorbed by ground [W/m2]
                                   tmp_rssun             , & ! out   - sunlit leaf stomatal resistance [s/m]
                                   tmp_rssha             , & ! out   - shaded leaf stomatal resistance [s/m]
                                   tmp_bgap              , & ! out   - between gap fraction [-]
                                   tmp_wgap              , & ! out   - within gap fraction [-]
                                   tmp_tgb               , & ! out   - bare ground temperature [K]
                                   tmp_tgv               , & ! out   - under canopy ground temperature [K]
                                   tmp_chv               , & ! out   - sensible heat exchange coefficient vegetated [-]
                                   tmp_chb               , & ! out   - sensible heat exchange coefficient bare-ground [-]
                                   tmp_shg               , & ! out   - veg ground sensible heat [+ to atm] [W/m2]
                                   tmp_shc               , & ! out   - canopy sensible heat [+ to atm] [W/m2]
                                   tmp_shb               , & ! out   - bare sensible heat [+ to atm] [W/m2]
                                   tmp_evg               , & ! out   - veg ground evaporation [+ to atm] [W/m2]
                                   tmp_evb               , & ! out   - bare soil evaporation [+ to atm] [W/m2]
                                   tmp_ghv               , & ! out   - veg ground heat flux [+ to soil] [W/m2]
                                   tmp_ghb               , & ! out   - bare ground heat flux [+ to soil] [W/m2]
                                   tmp_irg               , & ! out   - veg ground net LW radiation [+ to atm] [W/m2]
                                   tmp_irc               , & ! out   - canopy net LW radiation [+ to atm] [W/m2]
                                   tmp_irb               , & ! out   - bare net LW radiation [+ to atm] [W/m2]
                                   tmp_tr                , & ! out   - transpiration [ to atm] [W/m2]
                                   tmp_evc               , & ! out   - canopy evaporation heat [to atm] [W/m2]
                                   tmp_chleaf            , & ! out   - leaf exchange coefficient [-]
                                   tmp_chuc              , & ! out   - under canopy exchange coefficient [-]
                                   tmp_chv2              , & ! out   - veg 2m exchange coefficient [-]
                                   tmp_chb2              )   ! out   - bare 2m exchange coefficient [-]

            ! save state variables from local variables to global variables
            NOAHMP401_struc(n)%noahmp401(t)%sfcrunoff       = tmp_sfcrunoff
            NOAHMP401_struc(n)%noahmp401(t)%udrrunoff       = tmp_udrrunoff
            NOAHMP401_struc(n)%noahmp401(t)%smc(:)          = tmp_smc(:)
            NOAHMP401_struc(n)%noahmp401(t)%sh2o(:)         = tmp_sh2o(:)
            NOAHMP401_struc(n)%noahmp401(t)%tslb(:)         = tmp_tslb(:)
            NOAHMP401_struc(n)%noahmp401(t)%sneqv           = tmp_sneqv
            NOAHMP401_struc(n)%noahmp401(t)%snowh           = tmp_snowh
            NOAHMP401_struc(n)%noahmp401(t)%canwat          = tmp_canwat
            NOAHMP401_struc(n)%noahmp401(t)%acsnom          = tmp_acsnom
            NOAHMP401_struc(n)%noahmp401(t)%acsnow          = tmp_acsnow
            NOAHMP401_struc(n)%noahmp401(t)%isnow           = tmp_isnow
            NOAHMP401_struc(n)%noahmp401(t)%tv              = tmp_tv
            NOAHMP401_struc(n)%noahmp401(t)%tg              = tmp_tg
            NOAHMP401_struc(n)%noahmp401(t)%canice          = tmp_canice
            NOAHMP401_struc(n)%noahmp401(t)%canliq          = tmp_canliq
            NOAHMP401_struc(n)%noahmp401(t)%eah             = tmp_eah
            NOAHMP401_struc(n)%noahmp401(t)%tah             = tmp_tah
            NOAHMP401_struc(n)%noahmp401(t)%cm              = tmp_cm
            NOAHMP401_struc(n)%noahmp401(t)%ch              = tmp_ch
            NOAHMP401_struc(n)%noahmp401(t)%fwet            = tmp_fwet
            NOAHMP401_struc(n)%noahmp401(t)%sneqvo          = tmp_sneqvo
            NOAHMP401_struc(n)%noahmp401(t)%albold          = tmp_albold
            NOAHMP401_struc(n)%noahmp401(t)%qsnow           = tmp_qsnow
            NOAHMP401_struc(n)%noahmp401(t)%wslake          = tmp_wslake
            NOAHMP401_struc(n)%noahmp401(t)%zwt             = tmp_zwt
            NOAHMP401_struc(n)%noahmp401(t)%wa              = tmp_wa
            NOAHMP401_struc(n)%noahmp401(t)%wt              = tmp_wt
            NOAHMP401_struc(n)%noahmp401(t)%tsno(:)         = tmp_tsno(:)
            NOAHMP401_struc(n)%noahmp401(t)%zss(:)          = tmp_zss(:)
            NOAHMP401_struc(n)%noahmp401(t)%snowice(:)      = tmp_snowice(:)
            NOAHMP401_struc(n)%noahmp401(t)%snowliq(:)      = tmp_snowliq(:)
            NOAHMP401_struc(n)%noahmp401(t)%lfmass          = tmp_lfmass
            NOAHMP401_struc(n)%noahmp401(t)%rtmass          = tmp_rtmass
            NOAHMP401_struc(n)%noahmp401(t)%stmass          = tmp_stmass
            NOAHMP401_struc(n)%noahmp401(t)%wood            = tmp_wood
            NOAHMP401_struc(n)%noahmp401(t)%stblcp          = tmp_stblcp
            NOAHMP401_struc(n)%noahmp401(t)%fastcp          = tmp_fastcp
            NOAHMP401_struc(n)%noahmp401(t)%lai             = tmp_lai
            NOAHMP401_struc(n)%noahmp401(t)%sai             = tmp_sai
            NOAHMP401_struc(n)%noahmp401(t)%tauss           = tmp_tauss
            NOAHMP401_struc(n)%noahmp401(t)%smoiseq(:)      = tmp_smoiseq(:)
            NOAHMP401_struc(n)%noahmp401(t)%smcwtd          = tmp_smcwtd
            NOAHMP401_struc(n)%noahmp401(t)%deeprech        = tmp_deeprech
            NOAHMP401_struc(n)%noahmp401(t)%rech            = tmp_rech
            NOAHMP401_struc(n)%noahmp401(t)%grain           = tmp_grain
            NOAHMP401_struc(n)%noahmp401(t)%gdd             = tmp_gdd
            NOAHMP401_struc(n)%noahmp401(t)%pgs             = tmp_pgs
            NOAHMP401_struc(n)%noahmp401(t)%gecros_state(:) = tmp_gecros_state(:)

            ! save output variables from local variables to global variables
            NOAHMP401_struc(n)%noahmp401(t)%tsk       = tmp_tsk
            NOAHMP401_struc(n)%noahmp401(t)%hfx       = tmp_hfx
!           NOAHMP401_struc(n)%noahmp401(t)%fsh       = tmp_fsh

            NOAHMP401_struc(n)%noahmp401(t)%qfx       = tmp_qfx
            NOAHMP401_struc(n)%noahmp401(t)%lh        = tmp_lh
            NOAHMP401_struc(n)%noahmp401(t)%grdflx    = tmp_grdflx
            NOAHMP401_struc(n)%noahmp401(t)%albedo    = tmp_albedo
            NOAHMP401_struc(n)%noahmp401(t)%snowc     = tmp_snowc
            NOAHMP401_struc(n)%noahmp401(t)%emiss     = tmp_emiss
            NOAHMP401_struc(n)%noahmp401(t)%rs        = tmp_rs
            NOAHMP401_struc(n)%noahmp401(t)%t2mv      = tmp_t2mv
            NOAHMP401_struc(n)%noahmp401(t)%t2mb      = tmp_t2mb
            NOAHMP401_struc(n)%noahmp401(t)%q2mv      = tmp_q2mv
            NOAHMP401_struc(n)%noahmp401(t)%q2mb      = tmp_q2mb
            NOAHMP401_struc(n)%noahmp401(t)%trad      = tmp_trad
            NOAHMP401_struc(n)%noahmp401(t)%nee       = tmp_nee
            NOAHMP401_struc(n)%noahmp401(t)%gpp       = tmp_gpp
            NOAHMP401_struc(n)%noahmp401(t)%npp       = tmp_npp
            NOAHMP401_struc(n)%noahmp401(t)%fveg      = tmp_fveg
            NOAHMP401_struc(n)%noahmp401(t)%runsf     = tmp_runsf
            NOAHMP401_struc(n)%noahmp401(t)%runsb     = tmp_runsb
            NOAHMP401_struc(n)%noahmp401(t)%ecan      = tmp_ecan
! Direct soil evaporation does not include sublimation of the snowpack
! on the soil (by the strict ALMA definition of ESoil). - David Mocko
            NOAHMP401_struc(n)%noahmp401(t)%edir      = tmp_edir - tmp_subsnow
            NOAHMP401_struc(n)%noahmp401(t)%etran     = tmp_etran
            NOAHMP401_struc(n)%noahmp401(t)%rainf     = tmp_rainf
            NOAHMP401_struc(n)%noahmp401(t)%snowf     = tmp_snowf
            NOAHMP401_struc(n)%noahmp401(t)%fsa       = tmp_fsa
            NOAHMP401_struc(n)%noahmp401(t)%fira      = tmp_fira
            NOAHMP401_struc(n)%noahmp401(t)%apar      = tmp_apar
            NOAHMP401_struc(n)%noahmp401(t)%psn       = tmp_psn
            NOAHMP401_struc(n)%noahmp401(t)%sav       = tmp_sav
            NOAHMP401_struc(n)%noahmp401(t)%sag       = tmp_sag
            NOAHMP401_struc(n)%noahmp401(t)%rssun     = tmp_rssun
            NOAHMP401_struc(n)%noahmp401(t)%rssha     = tmp_rssha
            NOAHMP401_struc(n)%noahmp401(t)%bgap      = tmp_bgap
            NOAHMP401_struc(n)%noahmp401(t)%wgap      = tmp_wgap
            NOAHMP401_struc(n)%noahmp401(t)%tgb       = tmp_tgb
            NOAHMP401_struc(n)%noahmp401(t)%tgv       = tmp_tgv
            NOAHMP401_struc(n)%noahmp401(t)%chv       = tmp_chv
            NOAHMP401_struc(n)%noahmp401(t)%chb       = tmp_chb
            NOAHMP401_struc(n)%noahmp401(t)%shg       = tmp_shg
            NOAHMP401_struc(n)%noahmp401(t)%shc       = tmp_shc
            NOAHMP401_struc(n)%noahmp401(t)%shb       = tmp_shb
            NOAHMP401_struc(n)%noahmp401(t)%evg       = tmp_evg
            NOAHMP401_struc(n)%noahmp401(t)%evb       = tmp_evb
            NOAHMP401_struc(n)%noahmp401(t)%ghv       = tmp_ghv
            NOAHMP401_struc(n)%noahmp401(t)%ghb       = tmp_ghb
            NOAHMP401_struc(n)%noahmp401(t)%irg       = tmp_irg
            NOAHMP401_struc(n)%noahmp401(t)%irc       = tmp_irc
            NOAHMP401_struc(n)%noahmp401(t)%irb       = tmp_irb
            NOAHMP401_struc(n)%noahmp401(t)%tr        = tmp_tr
            NOAHMP401_struc(n)%noahmp401(t)%evc       = tmp_evc
            NOAHMP401_struc(n)%noahmp401(t)%chleaf    = tmp_chleaf
            NOAHMP401_struc(n)%noahmp401(t)%chuc      = tmp_chuc
            NOAHMP401_struc(n)%noahmp401(t)%chv2      = tmp_chv2
            NOAHMP401_struc(n)%noahmp401(t)%chb2      = tmp_chb2

            ![ 1] output variable: tsk (unit=K  ). ***  surface radiative temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RADT, value = NOAHMP401_struc(n)%noahmp401(t)%tsk, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 2] output variable: fsh (unit=W/m2). ***  sensible heat flux to atmosphere
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QH, value = NOAHMP401_struc(n)%noahmp401(t)%hfx, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

           !call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QH,value = NOAHMP401_struc(n)%noahmp401(t)%fsh,   &
           !                                  vlevel=1,unit="W m-2",direction="UP",surface_type=LIS_rc%lsm_index)


            ![ 3] output variable: lh (unit=W/m2). ***  latent heat flux
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QLE, value = NOAHMP401_struc(n)%noahmp401(t)%lh, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 4] output variable: grdflx (unit=W/m2). ***  ground/snow heat flux to soil
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QG, value = NOAHMP401_struc(n)%noahmp401(t)%grdflx, &
                                              vlevel=1, unit="W m-2", direction="DN", surface_type = LIS_rc%lsm_index)

            ![ 5] output variable: albedo (unit=- ). ***  total grid albedo
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO, value = NOAHMP401_struc(n)%noahmp401(t)%albedo, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 6] output variable: snowc (unit=-). ***  snow cover fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWCOVER, value = NOAHMP401_struc(n)%noahmp401(t)%snowc, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 7] output variable: smc (unit=m3/m3). ***  volumetric soil moisture
            do i=1, NOAHMP401_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = NOAHMP401_struc(n)%noahmp401(t)%smc(i), &
                                                  vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 8] output variable: sh2o (unit=m3/m3). ***  equilibrium volumetric liquid soil moisture content
            do i=1, NOAHMP401_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMLIQFRAC, value = NOAHMP401_struc(n)%noahmp401(t)%sh2o(i), &
                                                  vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 9] output variable: tslb (unit=K). ***  soil temperature
            do i=1, NOAHMP401_struc(n)%nsoil
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILTEMP, value = NOAHMP401_struc(n)%noahmp401(t)%tslb(i), &
                                                  vlevel=i, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 10] output variable: sneqv (unit=mm ). ***  snow water equivalent
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWE, value = NOAHMP401_struc(n)%noahmp401(t)%sneqv, &
                                              vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 11] output variable: snowh (unit=m ). ***  physical snow depth
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDEPTH, value = NOAHMP401_struc(n)%noahmp401(t)%snowh, &
                                              vlevel=1, unit="m ", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 12] output variable: canwat (unit=kg/m2). ***  total canopy water storage
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPINT, value = NOAHMP401_struc(n)%noahmp401(t)%canwat, &
                                              vlevel=1, unit="kg m-2", direction="- ", surface_type = LIS_rc%lsm_index)

            ![ 13] output variable: emiss (unit=- ). ***  surface bulk emissivity
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EMISSFORC, value = NOAHMP401_struc(n)%noahmp401(t)%emiss, &
                                              vlevel=1, unit="- ", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 14] output variable: rs (unit=s/m). ***  total stomatal resistance
            ! call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RS, value = NOAHMP401_struc(n)%noahmp401(t)%rs, &
            !                                  vlevel=1, unit="s/m", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 15] output variable: isnow (unit=-). ***  actual number of snow layers
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOWN_NLAYER, value = -1.0*NOAHMP401_struc(n)%noahmp401(t)%isnow, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 16] output variable: tv (unit=K ). ***  vegetation leaf temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_VEGT, value = NOAHMP401_struc(n)%noahmp401(t)%tv, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 17] output variable: tg (unit=K). ***  averaged ground surface temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GROUNDAVGT, value = NOAHMP401_struc(n)%noahmp401(t)%tg, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 18] output variable: canice (unit=mm). ***  canopy intercepted ice
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWEVEG, value = NOAHMP401_struc(n)%noahmp401(t)%canice, &
                                              vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 19] output variable: canliq (unit=mm). ***  canopy intercepted liquid water
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_INTL, value = NOAHMP401_struc(n)%noahmp401(t)%canliq, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 20] output variable: eah (unit=Pa  ). ***  canopy air vapor pressure
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_VP, value = NOAHMP401_struc(n)%noahmp401(t)%eah, &
                                              vlevel=1, unit="Pa", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 21] output variable: tah (unit=K  ). ***  canopy air temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_TEMP, value = NOAHMP401_struc(n)%noahmp401(t)%tah, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 22] output variable: cm (unit=s/m ). ***  bulk momentum drag coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CM, value = NOAHMP401_struc(n)%noahmp401(t)%cm, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 23] output variable: ch (unit=s/m ). ***  bulk sensible heat exchange coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CH, value = NOAHMP401_struc(n)%noahmp401(t)%ch, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 24] output variable: fwet (unit=-  ). ***  wetted or snowed fraction of canopy
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CANOPY_WF, value = NOAHMP401_struc(n)%noahmp401(t)%fwet, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 25] output variable: wslake (unit=mm). ***  lake water storage
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAKEWATER, value = NOAHMP401_struc(n)%noahmp401(t)%wslake, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 26] output variable: zwt (unit=m). ***  water table depth
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WATERTABLED, value = NOAHMP401_struc(n)%noahmp401(t)%zwt, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 27] output variable: wa (unit=mm). ***  water storage in the "aquifer"
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GWS, value = NOAHMP401_struc(n)%noahmp401(t)%wa, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 28] output variable: wt (unit=mm). ***  water in aquifer and saturated soil
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WT_AQUI_SATSOIL, value = NOAHMP401_struc(n)%noahmp401(t)%wt, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 29] output variable: tsno (unit=K). ***  snow layer temperature
            do i=1, NOAHMP401_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWTPROF, value = NOAHMP401_struc(n)%noahmp401(t)%tsno(i), &
                                                  vlevel=i, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 30] output variable: snowice (unit=mm ). ***  snow layer ice
            do i=1, NOAHMP401_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWICE, value = NOAHMP401_struc(n)%noahmp401(t)%snowice(i), &
                                                  vlevel=i, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 31] output variable: snowliq (unit=mm ). ***  snow layer liquid water
            do i=1, NOAHMP401_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWLIQ, value = NOAHMP401_struc(n)%noahmp401(t)%snowliq(i), &
                                                  vlevel=i, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ! Yeosang Yoon, for snow DA
            ! output variable: z_snow (unit=m). ***  snow layer-bottom depth from snow surface
            do i=1, NOAHMP401_struc(n)%nsnow
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW_LBDFSS, value = NOAHMP401_struc(n)%noahmp401(t)%zss(i), &
                                                  vlevel=i, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            end do

            ![ 32] output variable: lfmass (unit=g/m2). ***  leaf mass
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LEAFMASS, value = NOAHMP401_struc(n)%noahmp401(t)%lfmass, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 33] output variable: rtmass (unit=g/m2 ). ***  mass of fine roots
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ROOTMASS, value = NOAHMP401_struc(n)%noahmp401(t)%rtmass, &
                                              vlevel=1, unit="g m-2 ", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 34] output variable: stmass (unit=g/m2 ). ***  stem mass
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_STEMMASS, value = NOAHMP401_struc(n)%noahmp401(t)%stmass, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 35] output variable: wood (unit=g/m2). ***  mass of wood including woody roots
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WOODMASS, value = NOAHMP401_struc(n)%noahmp401(t)%wood, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 36] output variable: stblcp (unit=g/m2). ***  stable carbon in deep soil
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CARBON_DEEPSOIL, value = NOAHMP401_struc(n)%noahmp401(t)%stblcp, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 37] output variable: fastcp (unit=g/m2 ). ***  short-lived carbon in shallow soil
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CARBON_SHALLOWSOIL, value = NOAHMP401_struc(n)%noahmp401(t)%fastcp, &
                                              vlevel=1, unit="g m-2", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 38] output variable: lai (unit=-). ***  leave area index
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LAI, value = NOAHMP401_struc(n)%noahmp401(t)%lai, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 39] output variable: sai (unit=- ). ***  stem area index
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SAI, value = NOAHMP401_struc(n)%noahmp401(t)%sai, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 40] output variable: tauss (unit=- ). ***  snow aging factor
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWAGE, value = NOAHMP401_struc(n)%noahmp401(t)%tauss, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 41] output variable: smcwtd (unit=m3/m3). ***  soil moisture content in the layer to the water table when deep
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BETWEENWATER, value = NOAHMP401_struc(n)%noahmp401(t)%smcwtd, &
                                              vlevel=1, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 42] output variable: deeprech (unit=m). ***  recharge to the water table when groundwater is deep
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QRECTOGW, value = NOAHMP401_struc(n)%noahmp401(t)%deeprech, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 43] output variable: rech (unit=m). ***  recharge from the water table when groundwater is shallow
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QRECFROMGW, value = NOAHMP401_struc(n)%noahmp401(t)%rech, &
                                              vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 44] output variable: t2mv (unit=K). ***  2-m air temperature over vegetated part
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_VEGE2MT, value = NOAHMP401_struc(n)%noahmp401(t)%t2mv, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 45] output variable: t2mb (unit=K ). ***  2-m height air temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BARE2MT, value = NOAHMP401_struc(n)%noahmp401(t)%t2mb, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 46] output variable: q2mv (unit=kg/kg). ***  2-m mixing ratio of vegetation part
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_VEGE2MQ2, value = NOAHMP401_struc(n)%noahmp401(t)%q2mv, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 47] output variable: q2mb (unit=kg/kg). ***  2-m mixing ratio of bare ground part
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BARE2MQ2, value = NOAHMP401_struc(n)%noahmp401(t)%q2mb, &
                                              vlevel=1, unit="kg kg-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 48] output variable: trad (unit=K). ***  surface radiative temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RADT, value = NOAHMP401_struc(n)%noahmp401(t)%trad, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 49] output variable: nee (unit=g/m2/s ). ***  net ecosystem exchange of CO2
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_NEE, value = NOAHMP401_struc(n)%noahmp401(t)%nee, &
                                              vlevel=1, unit="g m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)

            ![ 50] output variable: gpp (unit=g/m2/s  ). ***  gross primary assimilation of carbon
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GPP, value = NOAHMP401_struc(n)%noahmp401(t)%gpp, &
                                              vlevel=1, unit="g m-2 s-1", direction="IN", surface_type = LIS_rc%lsm_index)

            ![ 51] output variable: npp (unit=g/m2/s). ***  net primary productivity of carbon
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_NPP, value = NOAHMP401_struc(n)%noahmp401(t)%npp, &
                                              vlevel=1, unit="g m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)

            ![ 52] output variable: fveg (unit=-). ***  Noah-MP green vegetation fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GREENNESS, value = NOAHMP401_struc(n)%noahmp401(t)%fveg, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 53] output variable: runsf (unit=mm/s). ***  surface runoff
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QS, value = NOAHMP401_struc(n)%noahmp401(t)%runsf, &
                                              vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)

            ![ 54] output variable: runsb (unit=mm/s ). ***  baseflow (saturation excess)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSB, value = NOAHMP401_struc(n)%noahmp401(t)%runsb, &
                                              vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)

            ![ 55] output variable: ecan (unit=mm/s ). ***  evaporation of intercepted water
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ECANOP, value = NOAHMP401_struc(n)%noahmp401(t)%ecan, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 56] output variable: edir (unit=mm/s ). ***  soil surface evaporation rate
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ESOIL, value = NOAHMP401_struc(n)%noahmp401(t)%edir, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 57] output variable: etran (unit=mm/s ). ***  transpiration rate
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_TVEG, value = NOAHMP401_struc(n)%noahmp401(t)%etran, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 58] output variable: fsa (unit=W/m2). ***  total absorbed solar radiation
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWNET, value = NOAHMP401_struc(n)%noahmp401(t)%fsa, &
                                              vlevel=1, unit="W m-2", direction="DN", surface_type = LIS_rc%lsm_index)

            ![ 59] output variable: fira (unit=W/m2 ). ***  total net longwave radiation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LWUP, value = NOAHMP401_struc(n)%noahmp401(t)%fira, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 60] output variable: apar (unit=W/m2). ***  photosynthesis active energy by canopy
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_APAR, value = NOAHMP401_struc(n)%noahmp401(t)%apar, &
                                              vlevel=1, unit="W m-2", direction="IN", surface_type = LIS_rc%lsm_index)

            ![ 61] output variable: psn (unit=umol/m2/s ). ***  total photosynthesis of CO2 [+]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PSCO2, value = NOAHMP401_struc(n)%noahmp401(t)%psn, &
                                              vlevel=1, unit="umol m-2 s-1", direction="IN", surface_type = LIS_rc%lsm_index)

            !![ 62] output variable: sav (unit=W/m2 ). ***  solar radiation absorbed by vegetation
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SAV, value = NOAHMP401_struc(n)%noahmp401(t)%sav, &
            !                                  vlevel=1, unit="W/m2 ", direction="IN", surface_type = LIS_rc%lsm_index)

            ![ 63] output variable: sag (unit=W/m2 ). ***  solar radiation absorbed by ground
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SAQ, value = NOAHMP401_struc(n)%noahmp401(t)%sag, &
            !                                  vlevel=1, unit="W/m2 ", direction="IN", surface_type = LIS_rc%lsm_index)

            ![ 64] output variable: rssun (unit=s/m). ***  sunlit leaf stomatal resistance
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RSSUN, value = NOAHMP401_struc(n)%noahmp401(t)%rssun, &
                                              vlevel=1, unit="s m-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 65] output variable: rssha (unit=s/m). ***  shaded leaf stomatal resistance
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RSSHA, value = NOAHMP401_struc(n)%noahmp401(t)%rssha, &
                                              vlevel=1, unit="s m-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 66] output variable: bgap (unit=-). ***  between gap fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BGAP, value = NOAHMP401_struc(n)%noahmp401(t)%bgap, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 67] output variable: wgap (unit=- ). ***  within gap fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WGAP, value = NOAHMP401_struc(n)%noahmp401(t)%wgap, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 68] output variable: tgb (unit=K). ***  bare ground temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_BARESOILT, value = NOAHMP401_struc(n)%noahmp401(t)%tgb, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 69] output variable: tgv (unit=K). ***  vegetated ground surface temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GROUNDVEGT, value = NOAHMP401_struc(n)%noahmp401(t)%tgv, &
                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 70] output variable: chv (unit=s/m). ***  sensible heat exchange coefficient over vegetated fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHV, value = NOAHMP401_struc(n)%noahmp401(t)%chv, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 71] output variable: chb (unit=s/m). ***  sensible heat exchange coefficient over bare-ground fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHB, value = NOAHMP401_struc(n)%noahmp401(t)%chb, &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 72] output variable: shg (unit=W/m2     ). ***  get ground sensible heat [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SHG, value = NOAHMP401_struc(n)%noahmp401(t)%shg, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 73] output variable: shc (unit=W/m2   ). ***  canopy sensible heat [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SHC, value = NOAHMP401_struc(n)%noahmp401(t)%shc, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 74] output variable: shb (unit=W/m2     ). ***  bare ground sensible heat [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SHB, value = NOAHMP401_struc(n)%noahmp401(t)%shb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 75] output variable: evg (unit=W/m2  ). ***  veg ground evaporation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVG, value = NOAHMP401_struc(n)%noahmp401(t)%evg, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 76] output variable: evb (unit=W/m2  ). ***  bare soil evaporation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVB, value = NOAHMP401_struc(n)%noahmp401(t)%evb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 77] output variable: ghv (unit=W/m2 ). ***  vegetated ground heat flux [+ to soil]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GHV, value = NOAHMP401_struc(n)%noahmp401(t)%ghv, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 78] output variable: ghb (unit=W/m2 ). ***  bare ground heat flux [+ to soil]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GHB, value = NOAHMP401_struc(n)%noahmp401(t)%ghb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 79] output variable: irg (unit=W/m2 ). ***  veg ground net long wave radiation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_IRV, value = NOAHMP401_struc(n)%noahmp401(t)%irg, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 80] output variable: irc (unit=W/m2 ). ***  canopy net long wave radiation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_IRC, value = NOAHMP401_struc(n)%noahmp401(t)%irc, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 81] output variable: irb (unit=W/m2 ). ***  bare net long wave radiation [+ to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_IRB, value = NOAHMP401_struc(n)%noahmp401(t)%irb, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 82] output variable: tr (unit=W/m2 ). ***  transpiration heat [to atm]
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_HTR, value = NOAHMP401_struc(n)%noahmp401(t)%tr, &
                                              vlevel=1, unit="W m-2", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 83] output variable: evc (unit=W/m2 ). ***  canopy evaporation heat [to atm]
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVC, value = NOAHMP401_struc(n)%noahmp401(t)%evc, &
            !                                  vlevel=1, unit="W/m2 ", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 84] output variable: chleaf (unit=m/s). ***  leaf exchange coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHLEAF, value = NOAHMP401_struc(n)%noahmp401(t)%chleaf, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 85] output variable: chuc (unit=m/s). ***  under canopy exchange coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHUC, value = NOAHMP401_struc(n)%noahmp401(t)%chuc, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 86] output variable: chv2 (unit=m/s). ***  veg 2m exchange coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHV2, value = NOAHMP401_struc(n)%noahmp401(t)%chv2, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 87] output variable: chb2 (unit=m/s). ***  bare 2m exchange coefficient
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CHB2, value = NOAHMP401_struc(n)%noahmp401(t)%chb2, &
                                              vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 88] output variable: evap (unit=kg/m2/s). ***  total evapotranspiration
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVAP, value = NOAHMP401_struc(n)%noahmp401(t)%qfx, &
                                              vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)

            ![ 89] output variable: rainf (unit=kg/m2). ***  precipitation rate
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RAINF, value = NOAHMP401_struc(n)%noahmp401(t)%rainf, &
                                              vlevel=1, unit="kg m-2 s-1", direction="DN", surface_type = LIS_rc%lsm_index)

            ![ 90] output variable: snowf (unit=kg/m2). ***  snowfall rate
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWF, value = NOAHMP401_struc(n)%noahmp401(t)%snowf, &
                                              vlevel=1, unit="kg m-2 s-1", direction="DN", surface_type = LIS_rc%lsm_index)

            ![ 91] LWnet
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_LWNET,vlevel=1,  &
                  value=(-1.0 * NOAHMP401_struc(n)%noahmp401(t)%fira), &
                  unit="W m-2", direction="DN", surface_type=LIS_rc%lsm_index)

            ! Code added by Zhuo Wang on 02/28/2019
            ![ 92] output variable: qsnbot (unit=kg m-2 s-1). ***  melting water out of snow bottom 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSM, value = tmp_qsnbot, &
                  vlevel=1, unit="kg m-2 s-1", direction="S2L", surface_type = LIS_rc%lsm_index)

            ![ 93] output variable: subsnow (unit=kg m-2 s-1). ***  snow sublimation 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SUBSNOW, value = tmp_subsnow, &
                  vlevel=1, unit="kg m-2 s-1", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 94] output variable: AvgSurfT (unit=K). *** average surface temperature 
            AvgSurfT_out = NOAHMP401_struc(n)%noahmp401(t)%fveg * NOAHMP401_struc(n)%noahmp401(t)%tv + &
                  (1.0-NOAHMP401_struc(n)%noahmp401(t)%fveg) * NOAHMP401_struc(n)%noahmp401(t)%tgb
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AVGSURFT, value = AvgSurfT_out, &
                  vlevel=1, unit="K", direction="-",surface_type = LIS_rc%lsm_index)

            ![ 95] TWS should be SWE + CanopInt + Soil moisture + WA - David Mocko
            TWS_out = NOAHMP401_struc(n)%noahmp401(t)%sneqv
            if ((NOAHMP401_struc(n)%noahmp401(t)%canliq.ge.0.0).and.   &
                (NOAHMP401_struc(n)%noahmp401(t)%canice.ge.0.0)) then
               TWS_out = TWS_out +                                     &
                            (NOAHMP401_struc(n)%noahmp401(t)%canliq  + &
                             NOAHMP401_struc(n)%noahmp401(t)%canice)
            endif
            do i = 1,NOAHMP401_struc(n)%nsoil
               TWS_out = TWS_out +                                     &
                      (NOAHMP401_struc(n)%noahmp401(t)%smc(i)  *       &
                      tmp_sldpth(i)*LIS_CONST_RHOFW)
            enddo
            if (NOAHMP401_struc(n)%noahmp401(t)%wa.ge.0.0) then
               TWS_out = TWS_out + NOAHMP401_struc(n)%noahmp401(t)%wa
            endif
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_TWS, value = TWS_out, &
                   vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

            ![ 96] Qa - Advective energy - Heat transferred to a snow cover by rain
            !         - (unit=W m-2) - added by David Mocko
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QA, value = tmp_pah, &
                  vlevel=1, unit="W m-2",direction="DN",surface_type=LIS_rc%lsm_index)

! Added water balance change terms - David Mocko
            endsm = 0.0
            do i = 1,tmp_nsoil
               endsm = endsm +                                         &
                          (tmp_smc(i) * tmp_sldpth(i) * LIS_CONST_RHOFW)
            enddo
            call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSOILMOIST,&
                     value=(endsm - startsm),vlevel=1,unit="kg m-2",   &
                     direction="INC",surface_type=LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSWE,      &
                     value=(NOAHMP401_struc(n)%noahmp401(t)%sneqv -    &
                            startswe),                                 &
                     vlevel=1,unit="kg m-2",direction="INC",           &
                     surface_type=LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELINTERCEPT,&
                     value=((NOAHMP401_struc(n)%noahmp401(t)%canliq +  &
                             NOAHMP401_struc(n)%noahmp401(t)%canice) - &
                            startint),                                 &
                     vlevel=1,unit="kg m-2",direction="INC",           &
                     surface_type=LIS_rc%lsm_index)
! For now, the ALMA standard does not provide a variable for the
! change in groundwater storage.  Instead, temporarily writing it
! to the DELSURFSTOR (which is the change in surface water storage).
! This is only a temporary fix, until LIS_MOC_DELGROUNDWATER or
! a similarly-named variable is added into LIS_histDataMod.F90.
            call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_DELSURFSTOR, &
                     value=(NOAHMP401_struc(n)%noahmp401(t)%wa -       &
                            startgw),                                  &
                     vlevel=1,unit="kg m-2",direction="INC",           &
                     surface_type=LIS_rc%lsm_index)

            ! reset forcing variables to zeros
            NOAHMP401_struc(n)%noahmp401(t)%tair = 0.0
            NOAHMP401_struc(n)%noahmp401(t)%psurf = 0.0
            NOAHMP401_struc(n)%noahmp401(t)%wind_e = 0.0
            NOAHMP401_struc(n)%noahmp401(t)%wind_n = 0.0
            NOAHMP401_struc(n)%noahmp401(t)%qair = 0.0
            NOAHMP401_struc(n)%noahmp401(t)%swdown = 0.0
            NOAHMP401_struc(n)%noahmp401(t)%lwdown = 0.0
            NOAHMP401_struc(n)%noahmp401(t)%prcp = 0.0

        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        NOAHMP401_struc(n)%forc_count = 0

    endif ! end of alarmCheck loop

    deallocate( tmp_sldpth )
    deallocate( tmp_shdfac_monthly )
    deallocate( tmp_soilcomp )
    deallocate( tmp_smc )
    deallocate( tmp_sh2o )
    deallocate( tmp_tslb )
    deallocate( tmp_tsno )
    deallocate( tmp_zss )
    deallocate( tmp_snowice )
    deallocate( tmp_snowliq )
    deallocate( tmp_smoiseq )
    deallocate( tmp_gecros_state )

end subroutine NoahMP401_main
