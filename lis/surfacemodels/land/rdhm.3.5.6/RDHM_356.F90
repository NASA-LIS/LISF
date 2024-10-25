!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine RDHM_356(Tair, Psurf, Wind_E, Wind_N, Qair,      &   ! Forcing
                    Rainf, Snowf, Swdown, Lwdown,           &   ! Forcing
                    Tair_min, Tair_max,                     &   ! Forcing/Parameter 
                    TempHeight, WindHeight,                 &   ! Observation heights of air temperature/humidity and wind
                    DT_SAC_SNOW17, DT_FRZ,                  &   ! Control Parameters: time steps
                    FRZ_VER_OPT, SACHTET_OPT,               &   ! Control parameters     
                    SNOW17_OPT, PET_OPT,                    &   ! Control parameters
                    NSTYP, NVTYP,                           &   ! Control parameter
                    PET_MON, PETADJ_MON, GRN_MON,           &   ! monthly PET climatology, monthly greenness
                    SoilAlb, SnowAlb,                       &   ! constant parameter
                    SOILTYP, VEGETYP,                       &   ! soil type and vegetation type 
                    UZTWM, UZFWM, UZK, PCTIM, ADIMP, RIVA,  &   ! SACPAR
                    ZPERC, REXP, LZTWM, LZFSM, LZFPM, LZSK, &   ! SACPAR
                    LZPK, PFREE, SIDE, RSERV, EFC,          &   ! SACPAR
                    TBOT, RSMAX, CKSL, ZBOT,                &   ! FRZPAR, STXT = SOILTYP * 1.0  
                    vegRCMIN, climRCMIN,                    &   ! vegetation lookup parameter  
                    RGL, HS, LAI, D50, CROOT, Z0, &   ! vegetation lookup parameter
                    CLAY, SAND, SATDK,                      &   ! soil lookup parameter          > these parameters in LIS configuration file
                    CZIL, FXEXP, vegRCMAX, TOPT, PC,        &   ! constant parameters           |
                    RDST, thresholdRCMIN,                   &   ! constant parameters          /
                    SFCREF, BAREADJ,                        &   ! constant parameters
                    UZTWC, UZFWC,LZTWC, LZFPC,LZFSC, ADIMC, &   ! SACHTET state variable (SATST)
                    TS0, TS1, TS2, TS3, TS4,                &   ! SACHTET, state variable (FRZST)
                    UZTWH, UZFWH,LZTWH, LZFSH, LZFPH,       &   ! SACHTET, state variable (FRZST)
                    SMC, SH2O,                              &   ! SACHTET state variable
                    NDINTW, NDSINT, DSINTW, DSINT,          &   ! SACHTET inputs 
                    NORMALIZE, ALON, ALAT,                  &   ! SACHTET inputs
                    SWINT, SWHINT, TSINT, FRZDUP, FRZDBT,   &   ! SACHTET output variables 
                    FROST, ALBEDO, SURF, GRND, TET, EDMND,  &   ! SACHTET output variables 
                    CH, CM,                                 &   ! SACHTET output varialbes
                    SCF, MFMAX, MFMIN, NMF, UADJ,           &   ! Snow-17 input parameters          
                    SI, MBASE, PXTEMP, PLWHC, TIPM, GM,     &   ! Snow-17 input parameters           
                    ELEV, LAEC, ADC, SNOW17_SWITCH,         &   ! Snow-17 input parameters 
                    WE, LIQW, NEGHS, TINDEX, ACCMAX,        &   ! Snow-17 state variables    
                    SNDPT, SNTMP, SB, SBAESC, SBWS,         &   ! Snow-17 state variables
                    STORAGE, AEADJ, EXLAG, NEXLAG,TA_PREV,  &   ! Snow-17 state variables 
                    SWE, SnowFrac, SnowDepth, RM)               ! Output of Snow-17
    
    use ESMF  
    use LIS_coreMod
    use LIS_FORC_AttributesMod
    use LIS_XMRG_Reader, only : latlon_to_hrap

    implicit none 

    real, intent(in)    :: Tair                ! air temperature [K]
    real, intent(in)    :: Psurf               ! surface air pressure [Pa]
    real, intent(in)    :: Wind_E              ! eastward wind [m s-1]
    real, intent(in)    :: Wind_N              ! northward wind [m s-1]
    real, intent(in)    :: Qair                ! near surface specific humidity [kg kg-1]
    real, intent(in)    :: Rainf               ! rainfall rate [kg m-2 s-1]
    real, intent(in)    :: Snowf               ! snowfall rate [kg m-2 s-1]
    real, intent(in)    :: Swdown              ! incident shortwave radiation [W m-2]
    real, intent(in)    :: Lwdown              ! incident longwave radiation [W m-2]
    real, intent(in)    :: Tair_min            ! spatial parameter, daily minimum air temperature [K]
    real, intent(in)    :: Tair_max            ! spatial parameter, daily maximum air temperature [K]
    real, intent(in)    :: TempHeight          ! observation height of temperature of humidity [m]
    real, intent(in)    :: WindHeight          ! observation height of wind [m]
    real, intent(in)    :: DT_SAC_SNOW17       ! simulation time interval of SAC model and Snow-17 [s]
    real, intent(in)    :: DT_FRZ              ! simulation time interval of frozen soil model [s]
    integer, intent(in) :: FRZ_VER_OPT         ! version number of frozen soil model. 1: old version, 2: new version
    integer, intent(in) :: SNOW17_OPT          ! option for snow-17. If SNOW17_OPT=1, use SNOW-17, otherwise, don't use
    integer, intent(in) :: SACHTET_OPT         ! option for SAC-HTET. If SACHTET_OPT=1, use SACHTET. otherwise, don't run 
    integer, intent(in) :: NSTYP               ! number of soil types
    integer, intent(in) :: NVTYP               ! number of vegetation types 
    integer, intent(in) :: NDINTW              ! number of desired soil layers for total and liquid soil moisture
    integer, intent(in) :: NDSINT              ! number of desired soil layers for soil temperature
    integer, intent(in) :: NORMALIZE           ! normalization flag for total and liquid soil moisture output (1-normalized, 0-not) 
    real, intent(in)    :: DSINTW(NDINTW)      ! thickness of desired soil layers for liquid and total soil moisture (cm)
    real, intent(in)    :: DSINT(NDSINT)       ! thickness of desired soil layers for soil temperature (cm)

    real, intent(in)    :: PET_MON(12)         ! multiband spatial parameter, monthly PET climatology, time series of spatial parameter, [mm] 
    real, intent(in)    :: GRN_MON(12)         ! multiband spatial parameter, monthly greenness climatology, time series of spatial parameter [-]
    real, intent(in)    :: PETADJ_MON(12)      ! lookup table parameter, adjustment of PET for 12 months 
    real, intent(in)    :: SoilAlb             ! snow free ALBEDO (default value 0.15)
    real, intent(in)    :: SnowAlb             ! snow ALBEDO (default value 0.7)

    integer, intent(in) :: SOILTYP
    integer, intent(in) :: VEGETYP
    ! sac parameters
    real, intent(in)    :: UZTWM               ! upper zone tension water maximum storage, mm
    real, intent(in)    :: UZFWM               ! upper zone free water maximum storage, mm
    real, intent(in)    :: UZK                 ! upper zone free water latent depletion rate, day^-1
    real, intent(in)    :: PCTIM               ! impervious fraction of the watershad area, -
    real, intent(in)    :: ADIMP               ! additional impervious area, -
    real, intent(in)    :: RIVA                ! riparian vegetation area, -
    real, intent(in)    :: ZPERC               ! maximum percolation rate, -
    real, intent(in)    :: REXP                ! 
    real, intent(in)    :: LZTWM               ! lower zone tension water maximum storage, mm
    real, intent(in)    :: LZFSM               ! lower zone supplemental free water (fast) maximum storage, mm
    real, intent(in)    :: LZFPM               ! lower zone primary free water (slow) maximum storage, mm
    real, intent(in)    :: LZSK                ! lower zone supplemental free water depletion rate, day^-1
    real, intent(in)    :: LZPK                ! lower zone primary free water depletion rate, day^-1
    real, intent(in)    :: PFREE               ! fraction percolation from upper to lower free water storage, -
    real, intent(in)    :: SIDE                ! ratio of deep recharge to channel base flow, -
    real, intent(in)    :: RSERV               ! fraction of lower zone free water not transferable to tension water,-
    real, intent(in)    :: EFC                 ! fraction of forest cover, -

    ! frozen soil parameters
    real, intent(in)    :: TBOT                ! spatial parameter, bottom boundary soil temperature, °C
    real, intent(in)    :: RSMAX               ! spatial parameter, maximum residual porosity (usually = 0.58)
    real, intent(in)    :: CKSL                ! spatial parameter, ratio of frozen to non-frozen surface (increase in frozen ground contact, usually = 8 s/m)
    real, intent(in)    :: ZBOT                ! spatial parameter, lower boundary depth (negative value, usually = -2.5 m)
    real, intent(in)    :: vegRCMIN            ! lookup table parameter, minimal stomatal resistance table for SACHTET, 14 values   
    real, intent(in)    :: climRCMIN           ! lookup table parameter, climate dependent miminal stomatal resistance for SACHTET, 14 values 
    real, intent(in)    :: RGL                 ! lookup table parameter, solar radiation threshold table for SACHTET, 14 values 
    real, intent(in)    :: HS                  ! lookup table parameter, vapor pressure resistance factor table for SACHTET, 14 values
    real, intent(in)    :: LAI                 ! lookup table parameter, leaf area index table for SACHTET, 14 values
    real, intent(in)    :: D50                 ! lookup table parameter, the depth (cm) table at which 50% roots are allocated for SACHTET, 14 values  
    real, intent(in)    :: CROOT               ! lookup table parameter, root distribution parameter table for SACHTET, 14 values
    real, intent(in)    :: Z0                  ! lookup table parameter, roughness coefficient of surface  
    real, intent(in)    :: CLAY                ! lookup table parameter, clay content for SACHTET, 12 values 
    real, intent(in)    :: SAND                ! lookup table parameter, sand content for sACHTET, 12 values
    real, intent(in)    :: SATDK               ! lookup table parameter, saturated hydraulic conductivityfor SACHTET, 12 values
    real, intent(in)    :: CZIL                ! constant parameter, default=0.12 Zilitinkevich 
    real, intent(in)    :: FXEXP               ! constant parameter, FXEXP(fxexp),(default=2.0) bare soil 
    real, intent(in)    :: vegRCMAX            ! constant parameter, RCMAX,(default=5000s/m) maximum stomatal resistance
    real, intent(in)    :: TOPT                ! constant parameter, TOPT,(default=298K)optimum air 
    real, intent(in)    :: PC                  ! constant parameter, plant coef. default pc = -1, 0.6 - 0.8
    integer, intent(in) :: PET_OPT             ! constant parameter, this constant allows selection of 
    integer, intent(in) :: RDST                ! constant parameter, default=1 means noah option,  this 
    real, intent(in)    :: thresholdRCMIN      ! constant parameter, this constant allows change of RCMIN (0.5) 
    real, intent(in)    :: SFCREF              ! constant parameter, reference wind speed for PET adjustment (4 m s-1)
    real, intent(in)    :: BAREADJ             ! constant parameter, Ek-Chen evaporation threshold switch. Bare soil evaporation option changes according to greenness. 
   
    !!!! Snow-17 input parameters 
    real, intent(in)    :: ALON                ! spatial parameter, logitude (-)
    real, intent(in)    :: ALAT                ! spatial parameter, latitude (-)
    real, intent(in)    :: SCF                 ! spatial parameter, snow fall correction factor (-)  
    real, intent(in)    :: MFMAX               ! spatial parameter, maximum melt factor (mm/(6hr°C))
    real, intent(in)    :: MFMIN               ! spatial parameter, minimum melt factor (mm/(6hr°C))
    real, intent(in)    :: NMF                 ! spatial parameter, maximum negative melt factor (mm/(6hr°C))
    real, intent(in)    :: UADJ                ! spatial parameter, the average wind function during rain-on-snow periods (mm mb-1) 
    real, intent(in)    :: SI                  ! spatial parameter, areal water-equivalent above which 100 percent areal snow cover 
    real, intent(in)    :: MBASE               ! spatial parameter, base temperature for non-rain melt factor (°C)
    real, intent(in)    :: PXTEMP              ! spatial parameter, temperature which spereates rain from snow (°C)
    real, intent(in)    :: PLWHC               ! spatial parameter, maximum amount of liquid-water held against gravity drainage 
    real, intent(in)    :: TIPM                ! spatial parameter, antecedent snow temperature index parameter (-)
    real, intent(in)    :: GM                  ! spatial parameter, daily ground melt (mm/day)
    real, intent(in)    :: ELEV                ! spatial parameter, elevation (m)
    real, intent(in)    :: LAEC                ! spatial parameter, snow-rain split temperature (°C)
    real, intent(in)    :: ADC(11)             ! multiband spatial parameter, Snow-17 curve coordinates
    integer, intent(in) :: SNOW17_SWITCH       ! constant parameter, switch variable change liquid water freezing version, 0: Victor's version, 1: Eric's version


    !!!! SACHTET state varialbes 
    real, intent(inout) :: UZTWC               ! upper zone tension water storage content, mm
    real, intent(inout) :: UZFWC               ! upper zone free water storage content, mm
    real, intent(inout) :: LZTWC               ! lower zone tension water storage content, mm
    real, intent(inout) :: LZFPC               ! lower zone primary free water storage content, mm
    real, intent(inout) :: LZFSC               ! lower zone supplemental free water storage content, mm
    real, intent(inout) :: ADIMC               ! additional impervious area content, mm
                                          
    real, intent(inout) :: TS0                 ! first soil layer temperature, °C
    real, intent(inout) :: TS1                 ! second soil layer temperature, °C 
    real, intent(inout) :: TS2                 ! third soil layer temperature, °C
    real, intent(inout) :: TS3                 ! fourth soil layer temperature, °C
    real, intent(inout) :: TS4                 ! fifth soil layer temperature, °C
    real, intent(inout) :: UZTWH               ! unfrozen upper zone tension water, mm
    real, intent(inout) :: UZFWH               ! unfrozen uppeer zone free water, mm
    real, intent(inout) :: LZTWH               ! unfrozen lower zone tension water, mm
    real, intent(inout) :: LZFSH               ! unfrozen lower zone supplemental free water, mm
    real, intent(inout) :: LZFPH               ! unfrozen lower zone primary free water, mm
                                          
    real, intent(inout) :: SMC(6)              ! volumetric content of total soil moisture at each layer   
    real, intent(inout) :: SH2O(6)             ! volumetric content of liquid soil moisture at each layer

    ! Snow-17 state variables
    real, intent(inout) :: WE                  ! snow water equivalent without liquid water (mm)
    real, intent(inout) :: LIQW                ! liquid water in snow (mm)
    real, intent(inout) :: NEGHS               ! negative snow heat (mm)
    real, intent(inout) :: TINDEX              ! antecedent temperature index (°C)
    real, intent(inout) :: ACCMAX              ! cumulated snow water including liquid (mm)
    real, intent(inout) :: SNDPT               ! snow depth (cm)
    real, intent(inout) :: SNTMP               ! average snow temperature (°C)
    real, intent(inout) :: SB                  ! the last highest snow water equivalent before any snow fall (mm)
    real, intent(inout) :: SBAESC              ! internal snow state during melt & new snow fall (cheched with Victor)
    real, intent(inout) :: SBWS                ! internal snow state during melt & new snow fall (checked with Victor)
    real, intent(inout) :: STORAGE             ! snow liquid water attanuation storage (mm)
    real, intent(inout) :: AEADJ               ! adjusted areal snow cover fraction (-)
    real, intent(inout) :: EXLAG(7)            ! array of lagged liquid water values (-) 
    integer, intent(inout) :: NEXLAG           ! number of ordinates in lagged liquid water array (EXLAG) (-)
    real, intent(inout) :: TA_PREV             ! air temperature of previous time step (°C)

    real, intent(out)   :: SWINT(NDINTW)       ! total volumetric soil moisture contents at desired soil layers (can be different from soil layers)
    real, intent(out)   :: SWHINT(NDINTW)      ! liquid volumetric soil moisture contents at desired soil layers (can be different from soil layers)
    real, intent(out)   :: TSINT(NDSINT)       ! soil temperature at desired soil layers (can be different from soil layers)
    real, intent(out)   :: FRZDUP              ! depth of the upper border of frozen ground from surface    m    <- cm
    real, intent(out)   :: FRZDBT              ! depth of the bottom border of frozen ground from surface   m    <- cm
    real, intent(out)   :: FROST               ! frost index -
    real, intent(out)   :: ALBEDO              ! land surface albedo 
    real, intent(out)   :: SURF                ! Qs      <=> SURF   simulated fast runoff (surface runoff?) mm s-1 <- mm/dt
    real, intent(out)   :: GRND                ! Qsb     <=> GRND   simulated slow runoff (baseflow?)       mm s-1 <- mm/dt
    real, intent(out)   :: TET                 ! Evap    <=> TET    simulated actual evapotranspiration     mm s-1 <- mm/dt
    real, intent(out)   :: EDMND               ! PotEvap <=> EDMND  potential evapotranspiration            mm s-1 <- mm
    real, intent(inout) :: CH                  ! surface layer exchage coefficient for heat and moisture [s/m]
    real, intent(inout) :: CM                  ! surface layer exchange coefficient for momentum (drag coefficient) [s/m]
   
    !!!! Snow-17 output variables
    real, intent(out)   :: SnowFrac            ! snow cover fraction, Snow-17 [-]
    real, intent(out)   :: SWE                 ! snow water equivalent, Snow-17 [kg m-2]
    real, intent(out)   :: SnowDepth           ! snow depth, Snow-17 [m]
    real, intent(out)   :: RM                  ! rain + melt [mm]
    
    ! local variables
    real :: TA                                 ! air temperature [C]    
    real :: PXV                                ! rainfall per one time interval [mm]
    real :: q2                                 ! specific humidity [kg kg-1]
    real :: q2sat                              ! saturated specific humidity
    real :: es                                 ! saturated air pressure [Pa]
    real :: solar                              ! downward shortwave radiation [W m-2]
    real :: rlwdn                              ! incoming long wave radiation [W m-2]
    real :: fdown                              ! downward total radiation (fdown=solar*(1.0-ALBEDO)+rlwdn) [W m-2]
    real :: sfcprs                             ! surface air pressure at reference height [Pa]
    real :: sfcspd, sfcref_l                  ! sufrface wind speed at reference height [m s-1]
    real :: WE_SAC                             ! snow water equivalent [mm]
    real :: AESC_SAC                           ! snow cover fraction [-]
    real :: SH_SAC                             ! snow depth in [cm]
    real :: dtday, DT                          ! time step of SAC model in day [day]
    real :: DTFRZ                              ! frozen soil model time step in seconds
    real :: ztmp                               ! air temperature observation height (default = 2m)
    real :: zspd                               ! wind speed observation height (defult = 2m)
    integer :: IDTFRZ                          ! number of time interval for ground model  
    integer :: YEAR, MONTH, DAY, HOUR, MINUTE  ! OHD time converted from LIS data structure
    real :: SACST_PRV(6)                       ! SAC-HTET state variables of previous time step
    real :: SACST(6)                           ! SAC-HTET state variables of current time step
    real :: SACPAR(20)                         ! SAC parameters
    real :: FRZST(10)                          ! frozen ground model state variables 
    real :: FRZPAR(15)                         ! frozen ground model parameters 
    !integer :: PENPT
    real :: PENPT
    integer :: IVERS
    
    real :: rtdis(5)                           ! root fraction in soil layers  
    real :: dksat                              ! saturated hydraulic conductivity 
    real :: dwsat                              ! saturated soil moisture diffusivity 
    real :: rcmin                              ! minimum stomatal resistance [s/m]
    real :: rcmax                              ! maximum stomatal resistance [s/m]
    real :: bexp 

    real :: rgl_t                              ! solar radiation threshold                            
    real :: hs_t                               ! vapor presure resistance factor 
    real :: xlai_t                             ! leaf area index 
    real :: z0_t                               ! roughness coefficient of surface 
    real :: sfsl                               ! daily sum of canopy resistance factors to solar radiation which will be used to distribute daily ET demand 
    real :: grn(12)                            ! monthly greenness fraction 
    real :: grn_ind                            ! climate index, average of monthly greenness fraction 
    real :: emiss 
    real :: pcinp
    real :: soil_alb, snow_alb 

    ! program calculated frozen soil parameters
    real :: STXT             !  -   soil texture (12 classes), -
    real :: SMAX             !  -   soil porosity / saturated soil moisture content
    real :: SFLD             !  -   field capacity / reference valumetric water content 
    real :: BRT              !  -   slope of the retension curve, CAMPBELL-COSBY soil retention curve power paramter
    real :: SRT              !  -   LOG(PSISAT)+BRT*LOG(SMAX)+2
    real :: QUARTZ           !  -   fraction of quartz in soil 
    real :: STYPE            !  -   soil type (11 - if coarser soil, 12 - if fine)
    integer:: NSOIL          !  -   total number of selected layers 
    integer:: NUPL           !  -   NUMBER OF UPPER ZONE LAYERS + RESIDULE LAYERS
    integer:: NSAC           !  -   number of soil layers to reflect SAC-SMA storage 
    integer:: NROOT          !  
    real :: RTUP             !  -   frzpar(6): upper layer depth adjustment, calculated by subroutine SOILPAR1 according to soiltype
    real :: RTLW             !  -   frzpar(7): lower layer depth adjustment, calculated by subroutine SOILPAR1 according to soiltype
    real :: PSISAT           ! kPa	frzpar(8): saturated soil pressure (from soil texture), calculated by subroutine SOILPAR1 according to soiltype 
    real :: SWLT             !	-	frzpar(9): wilting point (from soil texture), (prgram calculated), calculated by subroutine SOILPAR1 according to soiltype
    real :: ZSOIL(5)         !  m   first, secod, third, forth, and fifth layer depth as negative value 
    real :: SUPM             !  mm upper zone soil moisture storage, SOILPAR1 will convert the mm into m
    real :: SLWM             !  mm lower zone soil moisture storage, SOILPAR1 will convert the mm into m

    real :: d50x, crootx, rtcx, rtc, rdpth
    integer :: n
   
    real    :: CHANLOSS
    real    :: HRAPX  
    real    :: HRAPY  
    integer :: ERROR

    ! internal variables for Snow-17
    integer :: IDD           ! day
    integer :: IMN           ! month
    real    :: DT_snow17     ! simulation time step in hours
    real    :: TA_snow17     ! air temperature  °C
    real    :: PX            ! precipitation mm
    real    :: PCTS          ! snow fraction of precipitation
    real    :: TWE           ! total snow water, mm
    real    :: COVER         ! snow cover fraction
    real    :: DS            ! calculated using WE and SNDPT
    real    :: DTA           ! calculated using tairx and dtax
    real    :: DTAX
    real    :: PA            ! air pressure in mbar or hPa
    integer :: SWITCH        ! switch variable for liquid water freezing version: 0-Victor's verison, 1-Eric's version
    real    :: SNOF          ! the threshold of new snow fall after which model will leave snow depletion curve and assume a uniform snow cover
    integer :: dtlocal       ! difference between local time and Z time  
    integer :: dt_min        ! time step in mimutes
    real    :: pet_adj        

    real    :: solardt(144)  
    real    :: esat          ! external function 

    real :: offsetlocal  ! function to calculate dtlocal
    real :: rday  
    real :: MFMAXX, MFMINX, NMFX,  UADJX, TIPMX, GMX   ! recalculated 

    ! ESMF variables are used to adjust LIS time into OHD time 
    type(ESMF_Time) :: lis_time, ohd_time
    type(ESMF_TimeInterval) :: timeinterval1h
    integer :: rc 

    ! YEAR, MONTH, DAY and HOUR are used in solar radiation approximation if OHD forcing 
    ! data are used. OHD RDHM time is one hour behind LIS time. The following code does
    ! conversion from LIS time into RDHM time.  04/07/2014 Shugong Wang 

    ! initialize time 1 
    call ESMF_TimeSet(lis_time, yy = LIS_rc%yr, &
                                mm = LIS_rc%mo, &
                                dd = LIS_rc%da, &
                                h  = LIS_rc%hr, &
                                m  = LIS_rc%mn, s=0)
    
    ! set time interval of 1 hour
    call ESMF_TimeIntervalSet(timeinterval1h, h=1, rc=rc)
    ohd_time = lis_time - timeinterval1h
    call ESMF_TimeGet(ohd_time, yy=YEAR, mm=MONTH, dd=DAY, h=HOUR, m=MINUTE, rc=rc);
    rday = day * 1.0 

    
    ! calculate HRAPX and HRAPY 
    call latlon_to_hrap(ALON, ALAT, HRAPX, HRAPY) 
    HRAPX = floor(HRAPX+0.5)
    HRAPY = floor(HRAPY+0.5)
    ! snow-17 input
    if(SNOW17_OPT .eq. 1) then
        !!!! Snow-17 
        IDD = DAY
        IMN = MONTH
        DT_snow17 = DT_SAC_SNOW17/3600  ! from second to hour
        TA_snow17 = Tair - 273.16  
        !!! precipitation
        if (LIS_FORC_Snowf%selectOpt .eq. 1) then 
            PX  = (Rainf + Snowf) * DT_snow17 * 3600.0    ! from kg m-2 s-1 -> mm
        else
            PX  = Rainf * DT_snow17 * 3600.0              ! fron kg m-2 s-1 -> mm
        endif
        ! snow fraction of precpitation (PCTS). If both rain fall and snow fall are provided, 
        ! calculate PCTS. Otherwise, set PCTS = -1. Then PCTS is to be determined by Snow-17
        if ((LIS_FORC_Snowf%selectOpt .eq. 1) .and. (LIS_FORC_Rainf%selectOpt .eq. 1)) then
            PCTS = Snowf / (Rainf + Snowf) 
        else
            PCTS = -1.0
        endif
        DTAX = TA_PREV
        DTA = TA_snow17 - DTAX
        if(DTAX .gt. 0 .and. TA_snow17 .gt. 0) DTA = abs(DTA)
        if(DTAX .gt. 0 .and. TA_snow17 .lt. 0) DTA = TA_snow17
        if(DTAX .eq. TA_snow17) DTA = 0.0

        if(SNDPT .gt. 0.0) then
            DS = 0.1 * WE/SNDPT
        else
            DS = 0.1
        endif
       
        SNOF = 0.2 * DT_snow17 ! Zhengtao Cui: snof = 0.2f * dthr;
        
        if (LIS_FORC_Psurf%selectOpt .eq. 1) then
            PA = Psurf * 0.01 ! from Pa to hPa (mbar)
        else
            PA = (29.9 - 0.00335 * ELEV + 0.00022 *( 0.01 * ELEV )**2.4) * 33.86
        endif
        SWITCH = SNOW17_SWITCH
        
        ! convert unit from 6 hr to actual model time step for MFMAX, MFMIN
        MFMAXX = MFMAX / (6.0/DT_snow17)
        MFMINX = MFMIN / (6.0/DT_snow17)
        NMFX = NMF / (6.0/DT_snow17)
        UADJX = UADJ / (6.0/DT_snow17) 

        ! calculate TIPM again
        TIPMX = 1.0 - (1.0 - TIPM)**(DT_snow17/6.0)
        ! GM from mm/day into mm/dthr
        GMX = GM / (24.0/DT_snow17)

        ! call Snow-17 physics  
        call PACK19(IDD,        &    ! IDD = DAY   = LIS_rc%da 
                    IMN,        &    ! IMN = MONTH = LIS_rc%mo 
                    DT_snow17,  &    ! DT_SAC_SNOW17/3600 (from second to hour)
                    TA_snow17,  &    ! Tair - 273.16 (from K to °C)
                    PX,         &    ! PX = Rainf (+ Snowf) -> mm from kg m-2 s-1
                    PCTS,       &    ! PCTS = Snowf/(Rainf + Snowf) or -1.0
                    RM,         &    ! output variable
                    TWE,        &    ! output variable, TWE -> SWE
                    COVER,      &    ! output variable, Cover -> SnowFrac 
                    WE,         &    ! state variable
                    NEGHS,      &    ! state varialbe
                    LIQW,       &    ! state variable 
                    TINDEX,     &    ! state variable 
                    ACCMAX,     &    ! state variable 
                    SNDPT,      &    ! state variable 
                    SNTMP,      &    ! state variable 
                    SB,         &    ! state variable 
                    SBAESC,     &    ! state variable 
                    SBWS,       &    ! state variable
                    STORAGE,    &    ! state variable 
                    AEADJ,      &    ! state variable 
                    EXLAG,      &    ! state variable
                    NEXLAG,     &    ! state variable 
                    ALAT,       &    ! input parameter
                    SCF,        &    ! input parameter 
                    MFMAXX,     &    ! input parameter
                    MFMINX,     &    ! input parameter
                    NMFX,       &    ! input parameter
                    UADJX,      &    ! input parameter
                    SI,         &    ! input parameter
                    MBASE,      &    ! input parameter
                    PXTEMP,     &    ! input parameter
                    PLWHC,      &    ! input parameter 
                    TIPMX,      &    ! input parameter
                    GMX,        &    ! input parameter 
                    PA,         &    ! from Psurf or Elevation
                    LAEC,       &    ! input parameter
                    ADC,        &    ! input parameter
                    SNOF,       &    ! SNOF = 0.2 * DT_snow17
                    DS,         &    ! DS = 0.1 * WE/SNDPT or DS = 0.1 
                    DTA,        &    ! calculated
                    SWITCH)          ! SWITCH = SNOW17_SWITCH

        !Snow-17 uses TA_PREV in C  
        TA_PREV = Tair - 273.16  
        ! Snow-17 outputs:
        SWE       = TWE
        SnowFrac  = COVER
        SnowDepth = SNDPT * 0.01 ! from cm to m
    

        WE_SAC    = SWE                    ! kg m-2 <---> mm
        AESC_SAC  = SnowFrac               ! [-]
        SH_SAC    = SnowDepth * 100.0      ! from m to cm 
    else
        WE_SAC    = 0.0 
        AESC_SAC  = 0.0
        SH_SAC    = 0.0
    endif


    SUPM = UZTWM + UZFWM     ! upper zone tension water storage + upper zone free water storage
    SLWM = LZTWM + LZFSM + LZFPM ! lower zone tension water storage + lower zone supplemental free water storage + lower zone primary free water storage
    STXT = 1.0 * SOILTYP 
    call SOILPAR2(STXT, SAND, CLAY, SUPM, SLWM, SMAX,   & 
                  PSISAT, BRT, SWLT, QUARTZ, STYPE,     &
                  NSOIL, NUPL, NSAC, ZSOIL, RTUP, RTLW)    
   
    NROOT = NSOIL - 1

    ! calculate SFLD according to lines 729-738 of sachtetFuncBefore.cpp
    ! in sachtetFuncBefore, if(itxt<=2), itxt
    ! itxt = static_cast<int>(std::floor(stxt + 0.5f))-1) 
    ! So, in the fortran code, the threshold should be 3 instead 2 
    ! Shugong Wang, 03/24/2014 
    if(STXT .le. 3) then
        SFLD = SMAX * (10.0/PSISAT)**(-1.0/BRT)
    else
        SFLD = SMAX * (20.0/PSISAT)**(-1.0/BRT)
    endif

    ! calculate 
    ! dksat = SATDK(SOILTYP)
    dksat = SATDK
    ! dwsat = 0.102 * BRT * satdk(SOILTYP) * PSISAT / SMAX ! dwsat[pix]=0.102*brt[pix]*satdk[itxt]*psisat[pix]/smax[pix];
    dwsat = 0.102 * BRT * satdk * PSISAT / SMAX ! dwsat[pix]=0.102*brt[pix]*satdk[itxt]*psisat[pix]/smax[pix];
    bexp  = BRT                                 ! SLOPE OF THE RETENTION CURVE 
    ! time step processing
    dtday   = DT_SAC_SNOW17 / 86400.0       ! time step in Day, unit: s -> day
    DT      = dtday
    DTFRZ   = DT_FRZ                 ! unit: s -> s  
    IDTFRZ  = int(DT_SAC_SNOW17/DT_FRZ + 0.01)  ! number of interval of frozen model in SACHTET_time step
     

    ! assign parameters 
    ! plant coefficient
    pcinp = PC 
  
    ! monthly greenness fraction 
    grn(:) = GRN_MON(:)
    ! climate index 
    grn_ind = sum(grn)/12.0

    ! minimum stomatal resistance [s/m]
    if(grn_ind > 0.5) then
        !rcmin = climRCMIN(VEGETYP)
        rcmin = climRCMIN
    else
        !rcmin = vegRCMIN(VEGETYP)
        rcmin = vegRCMIN
    endif

    ! maximum stomatal resistance [s/m]
    rcmax = vegRCMAX

    ! solar radiation threshold
    !rgl_t  = RGL(VEGETYP)
    rgl_t  = RGL
    ! vapor pressure resistance factor
    !hs_t   = HS(VEGETYP)
    hs_t   = HS
    ! leaf are index
    !xlai_t = LAI(VEGETYP)    
    xlai_t = LAI
    ! roughness coefficient of land surface 
    !z0_t   = Z0(VEGETYP)
    z0_t   = Z0

    ! calculate root distribution, the code in sachtetFuncBefore.cpp should be checked with Victor and Zhengtao
    !crootx = CROOT(VEGETYP)
    crootx = CROOT
    !d50x   = D50(VEGETYP)
    d50x   = D50
   
    rtc = 0
    do n=1, NROOT-1
        ! temporary solution for zero d50x
        if(d50x .le. 0.001) d50x=0.001
        rtcx = 1.0 / (1.0 + (100.0*(ZSOIL(1)-ZSOIL(n+1))/d50x)**crootx) ! I need to check with Victor about the equaiton. 
        rtdis(n) = rtcx - rtc
        rtc = rtcx
    enddo
    rtdis(NROOT) = 1.0 - rtc

    
    CHANLOSS = 0
    !HRAPX    = 0  
    !HRAPY    = 0
    ERROR    = 0   

    !!! forcing data processing
    ! air temperature
    TA      = Tair - 273.16          ! from K to C, see sachet.c

    ! air pressure 
    if (LIS_FORC_Psurf%selectOpt .eq. 1) then
        sfcprs  = Psurf                  ! [Pa]
    else                                 
        ! calculated according to elevation from hPa/mbar to Pa
        sfcprs = (29.9 - 0.00335 * ELEV + 0.00022 *( 0.01 * ELEV )**2.4) * 33.86 * 100.0 
    endif
    
    ! wind speed 
    sfcspd  = sqrt(Wind_E*Wind_E + Wind_N*Wind_N)
   
    ! precipitation
    if(SNOW17_OPT .eq. 1) then
        ! if running SNOW-17 and SAC-HTET at the same time, replace
        ! precipitation from f2t with precipitation plus melt from
        ! SNOW-17, all contained in SNOW-17 RM variable
        PXV = RM
    else
        if (LIS_FORC_Snowf%selectOpt .eq. 1) then
            PXV = (Rainf + Snowf) * DT * 86400.0  ! from kg m-2 s-1 to mm/day
        else
            PXV = Rainf * DT * 86400.0            ! from kg m-2 s-1 to mm/day
        endif
    endif

    ! calculate the difference between local time and Z time
    !if(ALON .ge. 0.0) then
    !    dtlocal = int((ALON + 7.5)/15.0)
    !else
    !    dtlocal = -1 * int(abs(ALON - 7.5)/15.0)
    !endif
    ! it seens that dtlocal should be a postive number. Then the model should only be run 
    ! in the west hemisphere if using SOLTMXMN to approximate solar radiation. Shugong Wang
    ! 03/18/2014
    !if(ALON .le. 0.0) then
    !    dtlocal = floor((-ALON + 7.5)/15.0)
    !endif
    dtlocal = offsetlocal(HRAPX, HRAPY)

    ! shortwave radiation
    if (LIS_FORC_SWdown%selectOpt .eq. 1) then
        solar   = Swdown                 ! downward shortwave radiation [W m-2]
    else ! using OHD RDHM forcing, approximate solar radiation 
        call SOLTMXMN(MONTH, DAY*1.0, DT_SAC_SNOW17, dtlocal, Tair_max-273.16, Tair_min-273.16, &
                      ALAT, rcmin, rcmax, rgl_t, xlai_t, sfsl, solardt)
        solar = solardt(HOUR+1)          ! modification to match OHD RDHM time
    endif
    ! longwave radiation
    if (LIS_FORC_LWdown%selectOpt .eq. 1) then
        rlwdn = Lwdown
    else
        emiss = 1.0 - 0.261 * exp(-7.77E-4 * TA**2.0)       ! unit of TA is °C
        rlwdn = emiss * 5.672E-8 * Tair**4.0             ! unit of Tair is K
    endif
    
    if (SoilAlb .eq. -1.0) then 
        soil_alb = 0.15
    else
        soil_alb = SoilAlb
    endif

    if (SnowAlb .eq. -1.0) then
        snow_alb = 0.70
    else
        snow_alb = SoilAlb
    endif
    ALBEDO = soil_alb + SnowFrac*(snow_alb-soil_alb)
    fdown   = solar * (1.0 - ALBEDO) + rlwdn

    ! calculate q2 and q2sat
    es      = esat(Tair)             ! calculate saturated air pressure [Pa]
    q2sat   = 0.622 * es/(sfcprs-(1.0-0.622)*es)
    
    ! if Q2 is provided, overite q2 with Qair
    if(LIS_FORC_Qair%selectOpt .eq. 1) then
        q2 = Qair
    else
        ! this is the default setting in SACHTET, q2 is calculated. 
        call vapoprs(Tair, sfcprs, q2sat, q2) ! q2 is calculated in subroutine vapoprs
    endif
    ! set frozen model option
    IVERS = FRZ_VER_OPT
   

    ! setup observation heights for air temperature and wind
    ztmp  = TempHeight
    zspd  = WindHeight

    ! pack SAC state variables
    SACST(1) = UZTWC
    SACST(2) = UZFWC
    SACST(3) = LZTWC
    SACST(4) = LZFSC
    SACST(5) = LZFPC
    SACST(6) = ADIMC

    SACST_PRV = SACST

    ! pack SAC parameters
    SACPAR(1)  = UZTWM
    SACPAR(2)  = UZFWM
    SACPAR(3)  = UZK
    SACPAR(4)  = PCTIM
    SACPAR(5)  = ADIMP
    SACPAR(6)  = RIVA
    SACPAR(7)  = ZPERC
    SACPAR(8)  = REXP
    SACPAR(9)  = LZTWM
    SACPAR(10) = LZFSM
    SACPAR(11) = LZFPM
    SACPAR(12) = LZSK
    SACPAR(13) = LZPK
    SACPAR(14) = PFREE
    SACPAR(15) = SIDE
    SACPAR(16) = RSERV * (LZFSM + LZFPM)
    SACPAR(17) = EFC


    ! pack frozen ground model state variables 
    FRZST(1)  = TS0
    FRZST(2)  = TS1
    FRZST(3)  = TS2
    FRZST(4)  = TS3
    FRZST(5)  = TS4
    FRZST(6)  = UZTWH
    FRZST(7)  = UZFWH
    FRZST(8)  = LZTWH
    FRZST(9)  = LZFSH
    FRZST(10) = LZFPH

    ! pack frozen ground model parameters
    FRZPAR(1)  = STXT
    FRZPAR(2)  = TBOT
    FRZPAR(3)  = RSMAX
    FRZPAR(4)  = CKSL
    FRZPAR(5)  = ZBOT 
    FRZPAR(6)  = RTUP
    FRZPAR(7)  = RTLW
    FRZPAR(8)  = PSISAT
    FRZPAR(9)  = SWLT
    FRZPAR(10) = ZSOIL(1)
    FRZPAR(11) = ZSOIL(2)
    FRZPAR(12) = ZSOIL(3)
    FRZPAR(13) = ZSOIL(4)
    FRZPAR(14) = ZSOIL(5)

    ! set PET option. 
    PENPT = PET_OPT
    
    ! PENPT = -1 : physical Penman PET
    ! PENPT = 0  : climatological MONTHly PET
    ! PENPT = 1  : emperical Penman PET
    ! PENPT = 2, 3, ... : future support for other PET data product
    if(PENPT .eq. 0) then
        ! calculat sfsl, solardt seems a dummy argument
        call SOLTMXMN(MONTH, DAY*1.0, DT_SAC_SNOW17, dtlocal, Tair_max-273.16, Tair_min-273.16, &
                      ALAT, rcmin, rcmax, rgl_t, xlai_t, sfsl, solardt) 

        ! time step in minutes
        dt_min = int(DT_SAC_SNOW17/60.0)
        
        ! PET adjustment coefficient 
        ! pet_adj = PETADJ_MON(MONTH)

        ! After the following call, EDMND is daily PET. It will be converted 
        ! into PET of the time step in the call of HTET_FLAND1)
        call RDHM_GET_PED(YEAR, MONTH, DAY, dt_min, PET_MON, PETADJ_MON, EDMND) 
    endif
    ! set output variavles to 0 in case they are used in SAC-HTET
    TET = 0.0
    SURF = 0.0
    GRND = 0.0
   
    ! wind speed adjustment according to land cover type
    sfcspd = sfcspd * (log(ztmp/Z0)/log(zspd/Z0))
    sfcref_l = sfcref * (log(ztmp/Z0)/log(zspd/Z0))
    
    ! fix negative wind speed in SAC-HTET. 04/22/2014
    if(sfcspd .lt. 0) sfcspd = 0.1*sfcspd
    if(sfcref_l .lt. 0) sfcref_l = 0.1*sfcref_l
    ! fix problem caused by very small snow depth
    ! if snow depth is less than 0.1 mm, set it to 0 
    if(SH_SAC<0.01) then
        SH_SAC = 0.0 
        WE_SAC = 0.0
        AESC_SAC = 0.0 
    endif
    !!! call SAC-HTET physics 
    call HTET_FLAND1(PXV,        & !    (1)  PXV      = Rainf * DT * 86400.0   ! from kg m-2 s-1 to mm/day
                     EDMND,      & !    (2)  treat as output variable
                     TA,         & !    (3)  TA       = Tair - 273.16
                     WE_SAC,     & !    (4)  WE_SAC   = SWE or WE_SAC = 0
                     AESC_SAC,   & !    (5)  AESC_SAC = SnowFrac 
                     SH_SAC,     & !    (6)  SH_SAC   = SnowDepth * 100.0  
                     DT,         & !    (7)  DT = dtday = DT_SAC_SNOW17 / 86400.0 
                     SACST,      & !    (8)  pack ...
                     FRZST,      & !    (9)  pack ...
                     SACPAR,     & !    (10) pack ...
                     FRZPAR,     & !    (11) pack ...
                     NSOIL,      & !    (12) call SOILPAR1 
                     NUPL,       & !    (13) call SOILPAR1
                     NSAC,       & !    (14) call SOILPAR1
                     IVERS,      & !    (15) IVERS = FRZ_VER_OPT
                     SURF,       & !    (16) output variable
                     GRND,       & !    (17) output variable
                     TET,        & !    (18) output variable
                     SMC,        & !    (19) state variable
                     SH2O,       & !    (20) state variable
                     SACST_PRV,  & !    (21) SACST_PRV = SACST 
                     DTFRZ,      & !    (22) DTFRZ   = DT_FRZ
                     IDTFRZ,     & !    (23) IDTFRZ  = int(DT_SAC_SNOW17/DT_FRZ + 0.01)
                     FRZDUP,     & !    (24) output variable 
                     FRZDBT,     & !    (25) output variable
                     FROST,      & !    (26) output variable
                     TSINT,      & !    (27) output variable 
                     SWINT,      & !    (28) output variable
                     SWHINT,     & !    (29) output variable
                     DSINT,      & !    (30) input variable, thickness of desired soil layers for soil temperature (cm)
                     NDSINT,     & !    (31) input variable, number of desired soil layers for soil temperature
                     DSINTW,     & !    (32) input variable, thickness of desired soil layers for liquid and total soil moisture (cm)
                     NDINTW,     & !    (33) input variable, number of desired soil layers for total and liquid soil moisture
                     NORMALIZE,  & !    (34) input variable, normalization flag for total and liquid soil moisture output (1-normalized, 0-not) 
                     YEAR,       & !    (35) YEAR   = LIS_rc%yr
                     MONTH,      & !    (36) MONTH  = LIS_rc%mo
                     DAY,        & !    (37) DAY    = LIS_rc%da
                     HOUR,       & !    (38) HOUR   = LIS_rc%hr
                     fdown,      & !    (39) fdown   = solar * (1.0 - ALBEDO) + rlwdn
                     q2,         & !    (40) call vapoprs(Tair, sfcprs, q2sat, q2)
                     q2sat,      & !    (41) call vapoprs(Tair, sfcprs, q2sat, q2)
                     solar,      & !    (42) solar   = Swdown
                     nroot,      & !    (43) NROOT = NSOIL - 1
                     pcinp,      & !    (44) input parameter, pcinp = PC
                     smax,       & !    (45) call SOILPAR1, soil porosity
                     sfld,       & !    (46) SFLD = SMAX * (10.0/PSISAT)**(-1.0/BRT), or SFLD = SMAX * (20.0/PSISAT)**(-1.0/BRT)  field capacity
                     fxexp,      & !    (47) input parameter, (default=2.0) bare soil evaporation exponentin, FXEXP
                     rtdis,      & !    (48) ****
                     dwsat,      & !    (49) dwsat = 0.102 * BRT * satdk * PSISAT / SMAX 
                     dksat,      & !    (50) dksat = SATDK(SOILTYP)
                     bexp,       & !    (51) bexp  = BRT
                     rdst,       & !    (52) input constant variable 
                     sfcprs,     & !    (53) sfcprs  = Psurf
                     sfcspd,     & !    (54) sfcspd  = sqrt(Wind_E*Wind_E + Wind_N*Wind_N)
                     sfcref_l,   & !    (55) input parameter, reference wind speed for PET adjustment (m s-1) default 4 
                     rcmin,      & !    (56) rcmin = vegRCMIN(VEGETYP)
                     rcmax,      & !    (57) rcmax = vegRCMAX 
                     topt,       & !    (58) input constant parameter, optimal air temperature for transpiration (default=298K)
                     rgl_t,      & !    (59) solar radiation threshold, rgl_t  = RGL(VEGETYP)
                     hs_t,       & !    (60) vapor pressure resistance factor, hs_t   = HS(VEGETYP) 
                     xlai_t,     & !    (61) leaf area index, xlai_t = LAI(VEGETYP)    
                     sfsl,       & !    (62) ****   daily sum of canopy resistance factors to solar radiation which will be used to distribute daily ET demand 
                     grn,        & !    (63) input parameter, greenness fraction [0.0 - 1.0],  grn(:) = GRN_MON(:)
                     grn_ind,    & !    (64) climage index: average of monthly greenness,  grn_ind = sum(grn)/12.0   
                     ztmp,       & !    (65) input parameter,  ztmp  = TempHeight
                     z0_t,       & !    (67) lookup parameter, z0_t   = Z0(VEGETYP)
                     QUARTZ,     & !    (68) fraction of quartz in soil,  call SOILPAR2(...)
                     bareadj,    & !    (69) input parameter, BAREADJ 
                     czil,       & !    (70) input parameter, CZIL
                     cm,         & !    (71) output variable, surface layer exchange coefficient for momentum [s/m], calculated in htet_fland1
                     ch,         & !    (72) output variable, surface layer exchange coefficient for heat and moisture [s/m], calculated in htet_fland1
                     penpt,      & !    (73) input control parameter
                     CHANLOSS,   & !    (75) NEED SOME MORE INFORMATION
                     HRAPX,      & !    (76) NOT NEEDED
                     HRAPY,      & !    (77) NOT NEEDED
                     ERROR       ) !    (78) NOT NEEDED

    ! unpack SAC state variables
    UZTWC = SACST(1)
    UZFWC = SACST(2)
    LZTWC = SACST(3)
    LZFSC = SACST(4)
    LZFPC = SACST(5)
    ADIMC = SACST(6)

    ! unpack frozen ground model state variables 
    TS0   = FRZST(1) 
    TS1   = FRZST(2) 
    TS2   = FRZST(3) 
    TS3   = FRZST(4) 
    TS4   = FRZST(5) 
    UZTWH = FRZST(6) 
    UZFWH = FRZST(7) 
    LZTWH = FRZST(8) 
    LZFSH = FRZST(9) 
    LZFPH = FRZST(10)
end subroutine RDHM_356
