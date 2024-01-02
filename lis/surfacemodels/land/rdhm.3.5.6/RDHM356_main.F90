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
! !ROUTINE: RDHM356_main
! \label{RDHM356_main}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   3/7/14: Shugong Wang; initial implementation for RDHM356 with LIS-7
!
! !INTERFACE:
subroutine RDHM356_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod 
    use RDHM356_lsmMod
   !use other modules
  
    implicit none
! !ARGUMENTS:
    integer, intent(in)  :: n
    integer              :: t
    integer              :: i
    real                 :: dt
    real                 :: lat, lon
    integer              :: row, col
    logical              :: alarmCheck

!
! !DESCRIPTION:
!  This is the entry point for calling the RDHM356 physics.
!  This routine calls the {\tt RDHM\_356 } routine that performs the
!  land surface computations, to solve for water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for RDHM356
    real                 :: tmp_Tair               ! air temperature [K]
    real                 :: tmp_Psurf              ! surface air pressure [Pa]
    real                 :: tmp_Wind_E             ! eastward wind [m s-1]
    real                 :: tmp_Wind_N             ! northward wind [m s-1]
    real                 :: tmp_Qair               ! near surface specific humidity [kg kg-1]
    real                 :: tmp_Rainf              ! rainfall rate [kg m-2 s-1]
    real                 :: tmp_Snowf              ! snowfall rate [kg m-2 s-1]
    real                 :: tmp_Swdown             ! incident shortwave radiation [W m-2]
    real                 :: tmp_Lwdown             ! incident longwave radiation [W m-2]
    real                 :: tmp_Tair_min           ! daily minimum air temperature [K]
    real                 :: tmp_Tair_max           ! daily maximum air temperature [K]
    real                 :: tmp_TempHeight         ! observation height of temperature of humidity [m]
    real                 :: tmp_WindHeight         ! observation height of wind [m]
    real                 :: tmp_DT_SAC_SNOW17      ! simulation time interval of SAC model and Snow-17 [s]
    real                 :: tmp_DT_FRZ             ! simulation time interval of frozen soil model [s]
    integer              :: tmp_FRZ_VER_OPT        ! version number of frozen soil model. 1: old version, 2: new version [-]
    integer              :: tmp_SACHTET_OPT        ! option for snow-17. If SACHTET_OPT=1, run SACHTET, otherwise, don't run [-]
    integer              :: tmp_SNOW17_OPT         ! option for snow-17. If SNOW17_OPT=1, run SNOW-17, otherwise, don't run [-]
    integer              :: tmp_PET_OPT            ! if PET_OPT = 0, use non Penmann-based ETP;
                                                   !if penpt > 0 empirical Penmann equation; if penpt < 0, use energy based Pennman [-]
    integer              :: tmp_NSTYP              ! number of soil types [-]
    integer              :: tmp_NVTYP              ! number of vegetation types [-]
    real, allocatable    :: tmp_PET_MON(:)         ! multilevel monthly PET climatology, time series of spatial parameter [mm]
    real, allocatable    :: tmp_PETADJ_MON(:)      ! adjustment of PET for 12 months [-]
    real, allocatable    :: tmp_GRN_MON(:)         ! multilevel monthly greenness climatology, time series of spatial parameter [-] [-]
    real                 :: tmp_SoilAlb            ! snow free ALBEDO (default value 0.15) [-]
    real                 :: tmp_SnowAlb            ! snow ALBEDO (default value 0.7) [-]
    integer              :: tmp_SOILTYP            ! Soil type [-]
    integer              :: tmp_VEGETYP            ! Vegetation type [-]
    real                 :: tmp_UZTWM              ! upper zone tension water maximum storage [mm]
    real                 :: tmp_UZFWM              ! upper zone free water maximum storage [mm]
    real                 :: tmp_UZK                ! upper zone free water latent depletion rate [day^-1]
    real                 :: tmp_PCTIM              ! impervious fraction of the watershad area [-]
    real                 :: tmp_ADIMP              ! additional impervious area [-]
    real                 :: tmp_RIVA               ! riparian vegetation area [-]
    real                 :: tmp_ZPERC              ! maximum percolation rate [-]
    real                 :: tmp_REXP               ! exponent of the percolation equation (percolation parameter) [-]
    real                 :: tmp_LZTWM              ! lower zone tension water maximum storage [mm]
    real                 :: tmp_LZFSM              ! lower zone supplemental free water (fast) maximum storage [mm]
    real                 :: tmp_LZFPM              ! lower zone primary free water (slow) maximum storage [mm]
    real                 :: tmp_LZSK               ! lower zone supplemental free water depletion rate [day^-1]
    real                 :: tmp_LZPK               ! lower zone primary free water depletion rate [day^-1]
    real                 :: tmp_PFREE              ! fraction percolation from upper to lower free water storage [day^-1]
    real                 :: tmp_SIDE               ! ratio of deep recharge to channel base flow [-]
    real                 :: tmp_RSERV              ! fraction of lower zone free water not transferable to tension water [-]
    real                 :: tmp_EFC                ! fraction of forest cover [-]
    real                 :: tmp_TBOT               ! bottom boundary soil temperature [¡C]
    real                 :: tmp_RSMAX              ! maximum residual porosity (usually = 0.58) [-]
    real                 :: tmp_CKSL               ! ratio of frozen to non-frozen surface 
                                                   ! (increase in frozen ground contact, usually = 8 s/m) [s/m]
    real                 :: tmp_ZBOT               ! lower boundary depth (negative value, usually = -2.5 m) [m]
    real                 :: tmp_vegRCMIN           ! minimal stomatal resistance table for SACHTET, 14 values [s/m]
    real                 :: tmp_climRCMIN          ! climate dependent miminal stomatal resistance for SACHTET, 14 values [s/m]
    real                 :: tmp_RGL                ! solar radiation threshold table for SACHTET, 14 values [W m-2]
    real                 :: tmp_HS                 ! vapor pressure resistance factor table for SACHTET, 14 values [-]
    real                 :: tmp_LAI                ! leaf area index table for SACHTET, 14 values [-]
    real                 :: tmp_D50                ! the depth (cm) table at which 50% roots are allocated for SACHTET, 14 values [cm]
    real                 :: tmp_CROOT              ! root distribution parameter table for SACHTET, 14 values [-]
    real                 :: tmp_Z0                 ! roughness coefficient of surface [m]
    real                 :: tmp_CLAY               ! clay content for SACHTET, 12 values [-]
    real                 :: tmp_SAND               ! sand content for sACHTET, 12 values [-]
    real                 :: tmp_SATDK              ! saturated hydraulic conductivityfor SACHTET, 12 values [m s-1]
    real                 :: tmp_CZIL               ! default=0.12 Zilitinkevich [-]
    real                 :: tmp_FXEXP              ! FXEXP(fxexp),(default=2.0) bare soil [-]
    real                 :: tmp_vegRCMAX           ! RCMAX,(default=5000s/m) maximum stomatal resistance [s/m]
    real                 :: tmp_TOPT               ! TOPT,(default=298K) optimum air [K]
    real                 :: tmp_PC                 ! plant coef. default pc = -1, 0.6 - 0.8 [-]
    integer              :: tmp_RDST               ! default=1 means noah option,this constant allows
                                                   ! selection of tension water redistribution option, 
                                                   ! if rdst = 0 (ohd), use OHD version of SRT subroutine this SRT uses reference gradient instead an actual. 
                                                   ! if rdst = 1 ( noah), use Noah version of SRT subroutine [-]
    real                 :: tmp_thresholdRCMIN     ! this constant allows change of RCMIN (0.5) [s/m]
    real                 :: tmp_SFCREF             ! reference wind speed for PET adjustment (4 m s-1) [m/s]
    real                 :: tmp_BAREADJ            ! Ek-Chen evaporation threshold switch. Bare soil evaporation option changes according to greenness. [-]
    real                 :: tmp_UZTWC              ! upper zone tension water storage content [mm]
    real                 :: tmp_UZFWC              ! upper zone free water storage content [mm]
    real                 :: tmp_LZTWC              ! lower zone tension water storage content [mm]
    real                 :: tmp_LZFPC              ! lower zone primary free water storage content [mm]
    real                 :: tmp_LZFSC              ! lower zone supplemental free water storage content [mm]
    real                 :: tmp_ADIMC              ! additional impervious area content [mm]
    real                 :: tmp_TS0                ! first soil layer temperature [¡C]
    real                 :: tmp_TS1                ! second soil layer temperature [¡C]
    real                 :: tmp_TS2                ! third soil layer temperature [¡C]
    real                 :: tmp_TS3                ! fourth soil layer temperature [¡C]
    real                 :: tmp_TS4                ! fifth soil layer temperature [¡C]
    real                 :: tmp_UZTWH              ! unfrozen upper zone tension water [mm]
    real                 :: tmp_UZFWH              ! unfrozen uppeer zone free water [mm]
    real                 :: tmp_LZTWH              ! unfrozen lower zone tension water [mm]
    real                 :: tmp_LZFSH              ! unfrozen lower zone supplemental free water [mm]
    real                 :: tmp_LZFPH              ! unfrozen lower zone primary free water [mm]
    real, allocatable    :: tmp_SMC(:)             ! volumetric content of total soil moisture at each layer [m^3 m-3]
    real, allocatable    :: tmp_SH2O(:)            ! volumetric content of liquid soil moisture at each layer [m^3 m-3]
    integer              :: tmp_NDINTW             ! number of desired soil layers for total and liquid soil moisture [-]
    integer              :: tmp_NDSINT             ! number of desired soil layers for soil temperature [-]
    real, allocatable    :: tmp_DSINTW(:)          ! thickness of desired soil layers for liquid and total soil moisture [cm]
    real, allocatable    :: tmp_DSINT(:)           ! thickness of desired soil layers for soil temperature [cm]
    integer              :: tmp_NORMALIZE          ! normalization flag for total and liquid soil moisture output (1-normalized, 0-not) [-]
    real                 :: tmp_ALON               ! logitude [-]
    real                 :: tmp_ALAT               ! latitude [-]
    real, allocatable    :: tmp_SWINT(:)           ! total volumetric soil moisture contents at desired soil layers (can be different from soil layers) [-]
    real, allocatable    :: tmp_SWHINT(:)          ! liquid volumetric soil moisture contents at desired soil layers (can be different from soil layers) [-]
    real, allocatable    :: tmp_TSINT(:)           ! soil temperature at desired soil layers (can be different from soil layers) [-]
    real                 :: tmp_FRZDUP             ! depth of the upper border of frozen ground from surface [m]
    real                 :: tmp_FRZDBT             ! depth of the bottom border of frozen ground from surface [m]
    real                 :: tmp_FROST              ! frost index [-]
    real                 :: tmp_ALBEDO             ! land surface albedo [-]
    real                 :: tmp_SURF               ! Qs <=> SURF simulated fast runoff (surface runoff) [mm s-1]
    real                 :: tmp_GRND               ! Qsb <=> GRND simulated slow runoff (baseflow) [mm s-1]
    real                 :: tmp_TET                ! Evap <=> TET simulated actual evapotranspiration [mm s-1]
    real                 :: tmp_EDMND              ! PotEvap <=> EDMND potential evapotranspiration [mm s-1]
    real                 :: tmp_CH                 ! surface layer exchage coefficient for heat and moisture [s/m]
    real                 :: tmp_CM                 ! surface layer exchange coefficient for momentum (drag coefficient) [s/m]
    real                 :: tmp_SCF                ! snow fall correction factor [-]
    real                 :: tmp_MFMAX              ! maximum melt factor [mm/(6hr¡C)]
    real                 :: tmp_MFMIN              ! minimum melt factor [mm/(6hr¡C)]
    real                 :: tmp_NMF                ! maximum negative melt factor [mm/(6hr¡C)]
    real                 :: tmp_UADJ               ! the average wind function during rain-on-snow periods [mm/mb]
    real                 :: tmp_SI                 ! areal water-equivalent above which 100 percent areal snow cover [mm]
    real                 :: tmp_MBASE              ! base temperature for non-rain melt factor [¡C]
    real                 :: tmp_PXTEMP             ! temperature which spereates rain from snow [¡C]
    real                 :: tmp_PLWHC              ! maximum amount of liquid-water held against gravity drainage [-]
    real                 :: tmp_TIPM               ! antecedent snow temperature index parameter [-]
    real                 :: tmp_GM                 ! daily ground melt [mm/day]
    real                 :: tmp_ELEV               ! elevation [m]
    real                 :: tmp_LAEC               ! snow-rain split temperature [¡C]
    real, allocatable    :: tmp_ADC(:)             ! multilevel Snow-17 curve coordinates [-]
    integer              :: tmp_SNOW17_SWITCH      ! switch variable change liquid water freezing version, 0: Victor's version, 1: Eric's version [-]
    real                 :: tmp_WE                 ! snow water equivalent without liquid water [mm]
    real                 :: tmp_LIQW               ! liquid water in snow [mm]
    real                 :: tmp_NEGHS              ! negative snow heat [mm]
    real                 :: tmp_TINDEX             ! antecedent temperature index [¡C]
    real                 :: tmp_ACCMAX             ! cumulated snow water including liquid [mm]
    real                 :: tmp_SNDPT              ! snow depth [cm]
    real                 :: tmp_SNTMP              ! average snow temperature [¡C]
    real                 :: tmp_SB                 ! the last highest snow water equivalent before any snow fall [¡C]
    real                 :: tmp_SBAESC             ! internal snow state during melt & new snow fall (checked with Victor) [-]
    real                 :: tmp_SBWS               ! internal snow state during melt & new snow fall (checked with Victor) [-]
    real                 :: tmp_STORAGE            ! snow liquid water attenuation storage [mm]
    real                 :: tmp_AEADJ              ! adjusted areal snow cover fraction [-]
    real, allocatable    :: tmp_EXLAG(:)           ! array of lagged liquid water values [-]
    integer              :: tmp_NEXLAG             ! number of ordinates in lagged liquid water array (EXLAG) [-]
    real                 :: tmp_TA_PREV            ! air temperature of previous time step [-]
    real                 :: tmp_SWE                ! snow water equivalent, Snow-17 [kg m-2]
    real                 :: tmp_SnowFrac           ! snow cover fraction, Snow-17 [-]
    real                 :: tmp_SnowDepth          ! snow depth, Snow-17 [m]
    real                 :: tmp_RM                 ! rain + melt [mm]
    real                 :: SoilTemp(6)            ! soil temperature 
    allocate( tmp_PET_MON( 12 ) )
    allocate( tmp_PETADJ_MON( 12 ) )
    allocate( tmp_GRN_MON( 12 ) )
    allocate( tmp_SMC( 6 ) )
    allocate( tmp_SH2O( 6 ) )
    allocate( tmp_DSINTW( RDHM356_struc(n)%NDINTW ) )
    allocate( tmp_DSINT( RDHM356_struc(n)%NDSINT ) )
    allocate( tmp_SWINT( RDHM356_struc(n)%NDINTW ) )
    allocate( tmp_SWHINT( RDHM356_struc(n)%NDINTW ) )
    allocate( tmp_TSINT( RDHM356_struc(n)%NDSINT ) )
    allocate( tmp_ADC( 11 ) )
    allocate( tmp_EXLAG( 7 ) )

    ! check RDHM356 alarm. If alarm is ring, run model. 
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "RDHM356 model alarm")
    if (alarmCheck) Then
        do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            dt = LIS_rc%ts
            row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
            lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

            ! retrieve forcing data from RDHM356_struc(n)%rdhm356(t) and assign to local variables
            ! Tair: air temperature
            tmp_Tair       = RDHM356_struc(n)%rdhm356(t)%Tair   / RDHM356_struc(n)%forc_count
 
            ! Psurf: surface air pressure
            if(LIS_Forc_Psurf%selectOpt .eq. 1) then 
                tmp_Psurf  = RDHM356_struc(n)%rdhm356(t)%Psurf  / RDHM356_struc(n)%forc_count
            endif
 
            ! Wind_E: eastward wind
            tmp_Wind_E     = RDHM356_struc(n)%rdhm356(t)%Wind_E / RDHM356_struc(n)%forc_count
 
            ! Wind_N: northward wind
            tmp_Wind_N     = RDHM356_struc(n)%rdhm356(t)%Wind_N / RDHM356_struc(n)%forc_count
 
            ! Qair: near surface specific humidity
            if(LIS_Forc_Qair%selectOpt .eq. 1) then 
                tmp_Qair   = RDHM356_struc(n)%rdhm356(t)%Qair   / RDHM356_struc(n)%forc_count
            endif
 
            ! Rainf: rainfall rate
            tmp_Rainf      = RDHM356_struc(n)%rdhm356(t)%Rainf  / RDHM356_struc(n)%forc_count
 
            ! Snowf: snowfall rate
            if(LIS_Forc_Snowf%selectOpt .eq. 1) then 
                tmp_Snowf  = RDHM356_struc(n)%rdhm356(t)%Snowf  / RDHM356_struc(n)%forc_count
            endif
 
            ! Swdown: incident shortwave radiation
            tmp_Swdown     = RDHM356_struc(n)%rdhm356(t)%Swdown / RDHM356_struc(n)%forc_count
 
            ! Lwdown: incident longwave radiation
            if(LIS_Forc_Lwdown%selectOpt .eq. 1) then 
                tmp_Lwdown = RDHM356_struc(n)%rdhm356(t)%Lwdown / RDHM356_struc(n)%forc_count
            endif
 
            ! 
            ! check validity of Tair
            if(tmp_Tair .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable Tair in RDHM356"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Psurf
            if((LIS_Forc_Psurf%selectOpt .eq. 1) .and. (tmp_Psurf .eq. LIS_rc%udef)) then
                write(LIS_logunit, *) "undefined value found for forcing variable Psurf in RDHM356"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Wind_E
            if(tmp_Wind_E .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable Wind_E in RDHM356"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Wind_N
            if(tmp_Wind_N .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable Wind_N in RDHM356"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Qair
            if((LIS_Forc_Qair%selectOpt .eq. 1) .and. (tmp_Qair .eq. LIS_rc%udef)) then
                write(LIS_logunit, *) "undefined value found for forcing variable Qair in RDHM356"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Rainf
            if(tmp_Rainf .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable Rainf in RDHM356"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Snowf
            if((LIS_Forc_Snowf%selectOpt .eq. 1) .and. (tmp_Snowf .eq. LIS_rc%udef)) then
                write(LIS_logunit, *) "undefined value found for forcing variable Snowf in RDHM356"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Swdown
            if(tmp_Swdown .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable Swdown in RDHM356"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Lwdown
            if((LIS_Forc_Lwdown%selectOpt .eq. 1) .and. (tmp_Lwdown .eq. LIS_rc%udef)) then
                write(LIS_logunit, *) "undefined value found for forcing variable Lwdown in RDHM356"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
 
            ! get parameters 
            tmp_Tair_min               = RDHM356_struc(n)%rdhm356(t)%Tair_min                        
            tmp_Tair_max               = RDHM356_struc(n)%rdhm356(t)%Tair_max                        
            tmp_TempHeight             = RDHM356_struc(n)%TempHeight                      
            tmp_WindHeight             = RDHM356_struc(n)%WindHeight                      
            tmp_DT_SAC_SNOW17          = RDHM356_struc(n)%DT_SAC_SNOW17                   
            tmp_DT_FRZ                 = RDHM356_struc(n)%DT_FRZ                          
            tmp_FRZ_VER_OPT            = RDHM356_struc(n)%FRZ_VER_OPT                     
            tmp_SACHTET_OPT            = RDHM356_struc(n)%SACHTET_OPT                     
            tmp_SNOW17_OPT             = RDHM356_struc(n)%SNOW17_OPT                      
            tmp_PET_OPT                = RDHM356_struc(n)%PET_OPT                         
            tmp_NSTYP                  = RDHM356_struc(n)%NSTYP                           
            tmp_NVTYP                  = RDHM356_struc(n)%NVTYP                           
            tmp_PETADJ_MON(:)          = RDHM356_struc(n)%PETADJ_MON(:)                      
            tmp_GRN_MON(:)             = RDHM356_struc(n)%rdhm356(t)%GRN_MON(:)
            tmp_PET_MON(:)             = RDHM356_struc(n)%rdhm356(t)%PET_MON(:)
            tmp_ADC(:)                 = RDHM356_struc(n)%rdhm356(t)%ADC(:)
            tmp_SoilAlb                = RDHM356_struc(n)%rdhm356(t)%SoilAlb                         
            tmp_SnowAlb                = RDHM356_struc(n)%rdhm356(t)%SnowAlb                         
            tmp_SOILTYP                = RDHM356_struc(n)%rdhm356(t)%SOILTYP                         
            tmp_VEGETYP                = RDHM356_struc(n)%rdhm356(t)%VEGETYP                         
            tmp_UZTWM                  = RDHM356_struc(n)%rdhm356(t)%UZTWM                           
            tmp_UZFWM                  = RDHM356_struc(n)%rdhm356(t)%UZFWM                           
            tmp_UZK                    = RDHM356_struc(n)%rdhm356(t)%UZK                             
            tmp_PCTIM                  = RDHM356_struc(n)%rdhm356(t)%PCTIM                           
            tmp_ADIMP                  = RDHM356_struc(n)%rdhm356(t)%ADIMP                           
            tmp_RIVA                   = RDHM356_struc(n)%rdhm356(t)%RIVA                            
            tmp_ZPERC                  = RDHM356_struc(n)%rdhm356(t)%ZPERC                           
            tmp_REXP                   = RDHM356_struc(n)%rdhm356(t)%REXP                            
            tmp_LZTWM                  = RDHM356_struc(n)%rdhm356(t)%LZTWM                           
            tmp_LZFSM                  = RDHM356_struc(n)%rdhm356(t)%LZFSM                           
            tmp_LZFPM                  = RDHM356_struc(n)%rdhm356(t)%LZFPM                           
            tmp_LZSK                   = RDHM356_struc(n)%rdhm356(t)%LZSK                            
            tmp_LZPK                   = RDHM356_struc(n)%rdhm356(t)%LZPK                            
            tmp_PFREE                  = RDHM356_struc(n)%rdhm356(t)%PFREE                           
            tmp_SIDE                   = RDHM356_struc(n)%rdhm356(t)%SIDE                            
            tmp_RSERV                  = RDHM356_struc(n)%rdhm356(t)%RSERV                           
            tmp_EFC                    = RDHM356_struc(n)%rdhm356(t)%EFC                             
            tmp_TBOT                   = RDHM356_struc(n)%rdhm356(t)%TBOT                            
            tmp_RSMAX                  = RDHM356_struc(n)%rdhm356(t)%RSMAX                           
            tmp_CKSL                   = RDHM356_struc(n)%rdhm356(t)%CKSL                            
            tmp_ZBOT                   = RDHM356_struc(n)%rdhm356(t)%ZBOT                            
            tmp_vegRCMIN               = RDHM356_struc(n)%rdhm356(t)%vegRCMIN                        
            tmp_climRCMIN              = RDHM356_struc(n)%rdhm356(t)%climRCMIN                       
            tmp_RGL                    = RDHM356_struc(n)%rdhm356(t)%RGL                             
            tmp_HS                     = RDHM356_struc(n)%rdhm356(t)%HS                              
            tmp_LAI                    = RDHM356_struc(n)%rdhm356(t)%LAI                             
            tmp_D50                    = RDHM356_struc(n)%rdhm356(t)%D50                             
            tmp_CROOT                  = RDHM356_struc(n)%rdhm356(t)%CROOT                           
            tmp_Z0                     = RDHM356_struc(n)%rdhm356(t)%Z0                              
            tmp_CLAY                   = RDHM356_struc(n)%rdhm356(t)%CLAY                            
            tmp_SAND                   = RDHM356_struc(n)%rdhm356(t)%SAND                            
            tmp_SATDK                  = RDHM356_struc(n)%rdhm356(t)%SATDK                           
            tmp_CZIL                   = RDHM356_struc(n)%CZIL                            
            tmp_FXEXP                  = RDHM356_struc(n)%FXEXP                           
            tmp_vegRCMAX               = RDHM356_struc(n)%vegRCMAX                        
            tmp_TOPT                   = RDHM356_struc(n)%TOPT                            
            tmp_PC                     = RDHM356_struc(n)%PC                              
            tmp_RDST                   = RDHM356_struc(n)%RDST                            
            tmp_thresholdRCMIN         = RDHM356_struc(n)%thresholdRCMIN                  
            tmp_SFCREF                 = RDHM356_struc(n)%SFCREF                          
            tmp_BAREADJ                = RDHM356_struc(n)%BAREADJ                         
            tmp_NDINTW                 = RDHM356_struc(n)%NDINTW                          
            tmp_NDSINT                 = RDHM356_struc(n)%NDSINT                          
            tmp_DSINTW(:)              = RDHM356_struc(n)%DSINTW(:)                          
            tmp_DSINT(:)               = RDHM356_struc(n)%DSINT(:)                           
            tmp_NORMALIZE              = RDHM356_struc(n)%NORMALIZE                       
            tmp_ALON                   = RDHM356_struc(n)%rdhm356(t)%ALON                            
            tmp_ALAT                   = RDHM356_struc(n)%rdhm356(t)%ALAT                            
            tmp_SCF                    = RDHM356_struc(n)%rdhm356(t)%SCF                             
            tmp_MFMAX                  = RDHM356_struc(n)%rdhm356(t)%MFMAX                           
            tmp_MFMIN                  = RDHM356_struc(n)%rdhm356(t)%MFMIN                           
            tmp_NMF                    = RDHM356_struc(n)%rdhm356(t)%NMF                             
            tmp_UADJ                   = RDHM356_struc(n)%rdhm356(t)%UADJ                            
            tmp_SI                     = RDHM356_struc(n)%rdhm356(t)%SI                              
            tmp_MBASE                  = RDHM356_struc(n)%rdhm356(t)%MBASE                           
            tmp_PXTEMP                 = RDHM356_struc(n)%rdhm356(t)%PXTEMP                          
            tmp_PLWHC                  = RDHM356_struc(n)%rdhm356(t)%PLWHC                           
            tmp_TIPM                   = RDHM356_struc(n)%rdhm356(t)%TIPM                            
            tmp_GM                     = RDHM356_struc(n)%rdhm356(t)%GM                              
            tmp_ELEV                   = RDHM356_struc(n)%rdhm356(t)%ELEV                            
            tmp_LAEC                   = RDHM356_struc(n)%rdhm356(t)%LAEC                            
            tmp_SNOW17_SWITCH          = RDHM356_struc(n)%SNOW17_SWITCH                   
 
            ! get state variables
            tmp_UZTWC      = RDHM356_struc(n)%rdhm356(t)%UZTWC  
            tmp_UZFWC      = RDHM356_struc(n)%rdhm356(t)%UZFWC  
            tmp_LZTWC      = RDHM356_struc(n)%rdhm356(t)%LZTWC  
            tmp_LZFPC      = RDHM356_struc(n)%rdhm356(t)%LZFPC  
            tmp_LZFSC      = RDHM356_struc(n)%rdhm356(t)%LZFSC  
            tmp_ADIMC      = RDHM356_struc(n)%rdhm356(t)%ADIMC  
            tmp_TS0        = RDHM356_struc(n)%rdhm356(t)%TS0    
            tmp_TS1        = RDHM356_struc(n)%rdhm356(t)%TS1    
            tmp_TS2        = RDHM356_struc(n)%rdhm356(t)%TS2    
            tmp_TS3        = RDHM356_struc(n)%rdhm356(t)%TS3    
            tmp_TS4        = RDHM356_struc(n)%rdhm356(t)%TS4    
            tmp_UZTWH      = RDHM356_struc(n)%rdhm356(t)%UZTWH  
            tmp_UZFWH      = RDHM356_struc(n)%rdhm356(t)%UZFWH  
            tmp_LZTWH      = RDHM356_struc(n)%rdhm356(t)%LZTWH  
            tmp_LZFSH      = RDHM356_struc(n)%rdhm356(t)%LZFSH  
            tmp_LZFPH      = RDHM356_struc(n)%rdhm356(t)%LZFPH  
            tmp_SMC(:)     = RDHM356_struc(n)%rdhm356(t)%SMC(:)    
            tmp_SH2O(:)    = RDHM356_struc(n)%rdhm356(t)%SH2O(:)   
            tmp_WE         = RDHM356_struc(n)%rdhm356(t)%WE     
            tmp_LIQW       = RDHM356_struc(n)%rdhm356(t)%LIQW   
            tmp_NEGHS      = RDHM356_struc(n)%rdhm356(t)%NEGHS  
            tmp_TINDEX     = RDHM356_struc(n)%rdhm356(t)%TINDEX 
            tmp_ACCMAX     = RDHM356_struc(n)%rdhm356(t)%ACCMAX 
            tmp_SNDPT      = RDHM356_struc(n)%rdhm356(t)%SNDPT  
            tmp_SNTMP      = RDHM356_struc(n)%rdhm356(t)%SNTMP  
            tmp_SB         = RDHM356_struc(n)%rdhm356(t)%SB     
            tmp_SBAESC     = RDHM356_struc(n)%rdhm356(t)%SBAESC 
            tmp_SBWS       = RDHM356_struc(n)%rdhm356(t)%SBWS   
            tmp_STORAGE    = RDHM356_struc(n)%rdhm356(t)%STORAGE
            tmp_AEADJ      = RDHM356_struc(n)%rdhm356(t)%AEADJ  
            tmp_EXLAG(:)   = RDHM356_struc(n)%rdhm356(t)%EXLAG(:)  
            tmp_NEXLAG     = RDHM356_struc(n)%rdhm356(t)%NEXLAG 
            tmp_TA_PREV    = RDHM356_struc(n)%rdhm356(t)%TA_PREV
            tmp_CH         = RDHM356_struc(n)%rdhm356(t)%CH    
            tmp_CM         = RDHM356_struc(n)%rdhm356(t)%CM    
            
            if(tmp_CH>1.0)  tmp_CH = 1.0
            if(tmp_CM>1.0)  tmp_CM = 1.0
            ! 
            if(LIS_rc%tscount(n) .eq. 1) then
                tmp_TA_PREV = tmp_Tair - 273.16
            endif
            ! call model physics 
            call RDHM_356(tmp_Tair              , & ! in    - air temperature [K]
                          tmp_Psurf             , & ! in    - surface air pressure [Pa]
                          tmp_Wind_E            , & ! in    - eastward wind [m s-1]
                          tmp_Wind_N            , & ! in    - northward wind [m s-1]
                          tmp_Qair              , & ! in    - near surface specific humidity [kg kg-1]
                          tmp_Rainf             , & ! in    - rainfall rate [kg m-2 s-1]
                          tmp_Snowf             , & ! in    - snowfall rate [kg m-2 s-1]
                          tmp_Swdown            , & ! in    - incident shortwave radiation [W m-2]
                          tmp_Lwdown            , & ! in    - incident longwave radiation [W m-2]
                          tmp_Tair_min          , & ! in    - daily minimum air temperature [K]
                          tmp_Tair_max          , & ! in    - daily maximum air temperature [K]
                          tmp_TempHeight        , & ! in    - observation height of temperature of humidity [m]
                          tmp_WindHeight        , & ! in    - observation height of wind [m]
                          tmp_DT_SAC_SNOW17     , & ! in    - simulation time interval of SAC model and Snow-17 [s]
                          tmp_DT_FRZ            , & ! in    - simulation time interval of frozen soil model [s]
                          tmp_FRZ_VER_OPT       , & ! in    - version number of frozen soil model. 1: old version, 2: new version [-]
                          tmp_SACHTET_OPT       , & ! in    - option for snow-17. If SACHTET_OPT=1, run SACHTET, otherwise, don't run [-]
                          tmp_SNOW17_OPT        , & ! in    - option for snow-17. If SNOW17_OPT=1, run SNOW-17, otherwise, don't run [-]
                          tmp_PET_OPT           , & ! in    - if PET_OPT = 0, use non Penmann-based ETP;if penpt > 0 empirical Penmann equation; if penpt < 0, use energy based Pennman [-]
                          tmp_NSTYP             , & ! in    - number of soil types [-]
                          tmp_NVTYP             , & ! in    - number of vegetation types [-]
                          tmp_PET_MON           , & ! in    - multilevel monthly PET climatology, time series of spatial parameter [mm]
                          tmp_PETADJ_MON        , & ! in    - adjustment of PET for 12 months [-]
                          tmp_GRN_MON           , & ! in    - multilevel monthly greenness climatology, time series of spatial parameter [-] [-]
                          tmp_SoilAlb           , & ! in    - snow free ALBEDO (default value 0.15) [-]
                          tmp_SnowAlb           , & ! in    - snow ALBEDO (default value 0.7) [-]
                          tmp_SOILTYP           , & ! in    - Soil type [-]
                          tmp_VEGETYP           , & ! in    - Vegetation type [-]
                          tmp_UZTWM             , & ! in    - upper zone tension water maximum storage [mm]
                          tmp_UZFWM             , & ! in    - upper zone free water maximum storage [mm]
                          tmp_UZK               , & ! in    - upper zone free water latent depletion rate [day^-1]
                          tmp_PCTIM             , & ! in    - impervious fraction of the watershad area [-]
                          tmp_ADIMP             , & ! in    - additional impervious area [-]
                          tmp_RIVA              , & ! in    - riparian vegetation area [-]
                          tmp_ZPERC             , & ! in    - maximum percolation rate [-]
                          tmp_REXP              , & ! in    - exponent of the percolation equation (percolation parameter) [-]
                          tmp_LZTWM             , & ! in    - lower zone tension water maximum storage [mm]
                          tmp_LZFSM             , & ! in    - lower zone supplemental free water (fast) maximum storage [mm]
                          tmp_LZFPM             , & ! in    - lower zone primary free water (slow) maximum storage [mm]
                          tmp_LZSK              , & ! in    - lower zone supplemental free water depletion rate [day^-1]
                          tmp_LZPK              , & ! in    - lower zone primary free water depletion rate [day^-1]
                          tmp_PFREE             , & ! in    - fraction percolation from upper to lower free water storage [day^-1]
                          tmp_SIDE              , & ! in    - ratio of deep recharge to channel base flow [-]
                          tmp_RSERV             , & ! in    - fraction of lower zone free water not transferable to tension water [-]
                          tmp_EFC               , & ! in    - fraction of forest cover [-]
                          tmp_TBOT              , & ! in    - bottom boundary soil temperature [¡C]
                          tmp_RSMAX             , & ! in    - maximum residual porosity (usually = 0.58) [-]
                          tmp_CKSL              , & ! in    - ratio of frozen to non-frozen surface (increase in frozen ground contact, usually = 8 s/m) [s/m]
                          tmp_ZBOT              , & ! in    - lower boundary depth (negative value, usually = -2.5 m) [m]
                          tmp_vegRCMIN          , & ! in    - minimal stomatal resistance table for SACHTET, 14 values [s/m]
                          tmp_climRCMIN         , & ! in    - climate dependent miminal stomatal resistance for SACHTET, 14 values [s/m]
                          tmp_RGL               , & ! in    - solar radiation threshold table for SACHTET, 14 values [W m-2]
                          tmp_HS                , & ! in    - vapor pressure resistance factor table for SACHTET, 14 values [-]
                          tmp_LAI               , & ! in    - leaf area index table for SACHTET, 14 values [-]
                          tmp_D50               , & ! in    - the depth (cm) table at which 50% roots are allocated for SACHTET, 14 values [cm]
                          tmp_CROOT             , & ! in    - root distribution parameter table for SACHTET, 14 values [-]
                          tmp_Z0                , & ! in    - roughness coefficient of surface [m]
                          tmp_CLAY              , & ! in    - clay content for SACHTET, 12 values [-]
                          tmp_SAND              , & ! in    - sand content for sACHTET, 12 values [-]
                          tmp_SATDK             , & ! in    - saturated hydraulic conductivityfor SACHTET, 12 values [m s-1]
                          tmp_CZIL              , & ! in    - default=0.12 Zilitinkevich [-]
                          tmp_FXEXP             , & ! in    - FXEXP(fxexp),(default=2.0) bare soil [-]
                          tmp_vegRCMAX          , & ! in    - RCMAX,(default=5000s/m) maximum stomatal resistance [s/m]
                          tmp_TOPT              , & ! in    - TOPT,(default=298K) optimum air [K]
                          tmp_PC                , & ! in    - plant coef. default pc = -1, 0.6 - 0.8 [-]
                          tmp_RDST              , & ! in    - default=1 means noah option,this constant allows selection of tension water redistribution option, if rdst = 0 (ohd), use OHD version of SRT subroutine this SRT uses reference gradient instead an actual. if rdst = 1 ( noah), use Noah version of SRT subroutine [-]
                          tmp_thresholdRCMIN    , & ! in    - this constant allows change of RCMIN (0.5) [s/m]
                          tmp_SFCREF            , & ! in    - reference wind speed for PET adjustment (4 m s-1) [m/s]
                          tmp_BAREADJ           , & ! in    - Ek-Chen evaporation threshold switch. Bare soil evaporation option changes according to greenness. [-]
                          tmp_UZTWC             , & ! inout - upper zone tension water storage content [mm]
                          tmp_UZFWC             , & ! inout - upper zone free water storage content [mm]
                          tmp_LZTWC             , & ! inout - lower zone tension water storage content [mm]
                          tmp_LZFPC             , & ! inout - lower zone primary free water storage content [mm]
                          tmp_LZFSC             , & ! inout - lower zone supplemental free water storage content [mm]
                          tmp_ADIMC             , & ! inout - additional impervious area content [mm]
                          tmp_TS0               , & ! inout - first soil layer temperature [¡C]
                          tmp_TS1               , & ! inout - second soil layer temperature [¡C]
                          tmp_TS2               , & ! inout - third soil layer temperature [¡C]
                          tmp_TS3               , & ! inout - fourth soil layer temperature [¡C]
                          tmp_TS4               , & ! inout - fifth soil layer temperature [¡C]
                          tmp_UZTWH             , & ! inout - unfrozen upper zone tension water [mm]
                          tmp_UZFWH             , & ! inout - unfrozen uppeer zone free water [mm]
                          tmp_LZTWH             , & ! inout - unfrozen lower zone tension water [mm]
                          tmp_LZFSH             , & ! inout - unfrozen lower zone supplemental free water [mm]
                          tmp_LZFPH             , & ! inout - unfrozen lower zone primary free water [mm]
                          tmp_SMC               , & ! inout - volumetric content of total soil moisture at each layer [m^3 m-3]
                          tmp_SH2O              , & ! inout - volumetric content of liquid soil moisture at each layer [m^3 m-3]
                          tmp_NDINTW            , & ! in    - number of desired soil layers for total and liquid soil moisture [-]
                          tmp_NDSINT            , & ! in    - number of desired soil layers for soil temperature [-]
                          tmp_DSINTW            , & ! in    - thickness of desired soil layers for liquid and total soil moisture [cm]
                          tmp_DSINT             , & ! in    - thickness of desired soil layers for soil temperature [cm]
                          tmp_NORMALIZE         , & ! in    - normalization flag for total and liquid soil moisture output (1-normalized, 0-not) [-]
                          tmp_ALON              , & ! in    - logitude [-]
                          tmp_ALAT              , & ! in    - latitude [-]
                          tmp_SWINT             , & ! out   - total volumetric soil moisture contents at desired soil layers (can be different from soil layers) [-]
                          tmp_SWHINT            , & ! out   - liquid volumetric soil moisture contents at desired soil layers (can be different from soil layers) [-]
                          tmp_TSINT             , & ! out   - soil temperature at desired soil layers (can be different from soil layers) [-]
                          tmp_FRZDUP            , & ! out   - depth of the upper border of frozen ground from surface [m]
                          tmp_FRZDBT            , & ! out   - depth of the bottom border of frozen ground from surface [m]
                          tmp_FROST             , & ! out   - frost index [-]
                          tmp_ALBEDO            , & ! out   - land surface albedo [-]
                          tmp_SURF              , & ! out   - Qs <=> SURF simulated fast runoff (surface runoff) [mm s-1]
                          tmp_GRND              , & ! out   - Qsb <=> GRND simulated slow runoff (baseflow) [mm s-1]
                          tmp_TET               , & ! out   - Evap <=> TET simulated actual evapotranspiration [mm s-1]
                          tmp_EDMND             , & ! out   - PotEvap <=> EDMND potential evapotranspiration [mm s-1]
                          tmp_CH                , & ! inout - surface layer exchage coefficient for heat and moisture [s/m]
                          tmp_CM                , & ! inout  - surface layer exchange coefficient for momentum (drag coefficient) [s/m]
                          tmp_SCF               , & ! in    - snow fall correction factor [-]
                          tmp_MFMAX             , & ! in    - maximum melt factor [mm/(6hr¡C)]
                          tmp_MFMIN             , & ! in    - minimum melt factor [mm/(6hr¡C)]
                          tmp_NMF               , & ! in    - maximum negative melt factor [mm/(6hr¡C)]
                          tmp_UADJ              , & ! in    - the average wind function during rain-on-snow periods [mm/mb]
                          tmp_SI                , & ! in    - areal water-equivalent above which 100 percent areal snow cover [mm]
                          tmp_MBASE             , & ! in    - base temperature for non-rain melt factor [¡C]
                          tmp_PXTEMP            , & ! in    - temperature which spereates rain from snow [¡C]
                          tmp_PLWHC             , & ! in    - maximum amount of liquid-water held against gravity drainage [-]
                          tmp_TIPM              , & ! in    - antecedent snow temperature index parameter [-]
                          tmp_GM                , & ! in    - daily ground melt [mm/day]
                          tmp_ELEV              , & ! in    - elevation [m]
                          tmp_LAEC              , & ! in    - snow-rain split temperature [¡C]
                          tmp_ADC               , & ! in    - multilevel Snow-17 curve coordinates [-]
                          tmp_SNOW17_SWITCH     , & ! in    - switch variable change liquid water freezing version, 0: Victor's version, 1: Eric's version [-]
                          tmp_WE                , & ! inout - snow water equivalent without liquid water [mm]
                          tmp_LIQW              , & ! inout - liquid water in snow [mm]
                          tmp_NEGHS             , & ! inout - negative snow heat [mm]
                          tmp_TINDEX            , & ! inout - antecedent temperature index [¡C]
                          tmp_ACCMAX            , & ! inout - cumulated snow water including liquid [mm]
                          tmp_SNDPT             , & ! inout - snow depth [cm]
                          tmp_SNTMP             , & ! inout - average snow temperature [¡C]
                          tmp_SB                , & ! inout - the last highest snow water equivalent before any snow fall [¡C]
                          tmp_SBAESC            , & ! inout - internal snow state during melt & new snow fall (checked with Victor) [-]
                          tmp_SBWS              , & ! inout - internal snow state during melt & new snow fall (checked with Victor) [-]
                          tmp_STORAGE           , & ! inout - snow liquid water attenuation storage [mm]
                          tmp_AEADJ             , & ! inout - adjusted areal snow cover fraction [-]
                          tmp_EXLAG             , & ! inout - array of lagged liquid water values [-]
                          tmp_NEXLAG            , & ! inout - number of ordinates in lagged liquid water array (EXLAG) [-]
                          tmp_TA_PREV           , & ! inout - air temperature of previous time step [-]
                          tmp_SWE               , & ! out   - snow water equivalent, Snow-17 [kg m-2]
                          tmp_SnowFrac          , & ! out   - snow cover fraction, Snow-17 [-]
                          tmp_SnowDepth         , & ! out   - snow depth, Snow-17 [m]
                          tmp_RM                )   ! out   - rain + melt [mm]
 
    
            ! save state variables from local variables to global variables
            RDHM356_struc(n)%rdhm356(t)%UZTWC      = tmp_UZTWC     
            RDHM356_struc(n)%rdhm356(t)%UZFWC      = tmp_UZFWC     
            RDHM356_struc(n)%rdhm356(t)%LZTWC      = tmp_LZTWC     
            RDHM356_struc(n)%rdhm356(t)%LZFPC      = tmp_LZFPC     
            RDHM356_struc(n)%rdhm356(t)%LZFSC      = tmp_LZFSC     
            RDHM356_struc(n)%rdhm356(t)%ADIMC      = tmp_ADIMC     
            RDHM356_struc(n)%rdhm356(t)%TS0        = tmp_TS0       
            RDHM356_struc(n)%rdhm356(t)%TS1        = tmp_TS1       
            RDHM356_struc(n)%rdhm356(t)%TS2        = tmp_TS2       
            RDHM356_struc(n)%rdhm356(t)%TS3        = tmp_TS3       
            RDHM356_struc(n)%rdhm356(t)%TS4        = tmp_TS4       
            RDHM356_struc(n)%rdhm356(t)%UZTWH      = tmp_UZTWH     
            RDHM356_struc(n)%rdhm356(t)%UZFWH      = tmp_UZFWH     
            RDHM356_struc(n)%rdhm356(t)%LZTWH      = tmp_LZTWH     
            RDHM356_struc(n)%rdhm356(t)%LZFSH      = tmp_LZFSH     
            RDHM356_struc(n)%rdhm356(t)%LZFPH      = tmp_LZFPH     
            RDHM356_struc(n)%rdhm356(t)%SMC(:)     = tmp_SMC(:)    
            RDHM356_struc(n)%rdhm356(t)%SH2O(:)    = tmp_SH2O(:)   
            RDHM356_struc(n)%rdhm356(t)%WE         = tmp_WE        
            RDHM356_struc(n)%rdhm356(t)%LIQW       = tmp_LIQW      
            RDHM356_struc(n)%rdhm356(t)%NEGHS      = tmp_NEGHS     
            RDHM356_struc(n)%rdhm356(t)%TINDEX     = tmp_TINDEX    
            RDHM356_struc(n)%rdhm356(t)%ACCMAX     = tmp_ACCMAX    
            RDHM356_struc(n)%rdhm356(t)%SNDPT      = tmp_SNDPT     
            RDHM356_struc(n)%rdhm356(t)%SNTMP      = tmp_SNTMP     
            RDHM356_struc(n)%rdhm356(t)%SB         = tmp_SB        
            RDHM356_struc(n)%rdhm356(t)%SBAESC     = tmp_SBAESC    
            RDHM356_struc(n)%rdhm356(t)%SBWS       = tmp_SBWS      
            RDHM356_struc(n)%rdhm356(t)%STORAGE    = tmp_STORAGE   
            RDHM356_struc(n)%rdhm356(t)%AEADJ      = tmp_AEADJ     
            RDHM356_struc(n)%rdhm356(t)%EXLAG(:)   = tmp_EXLAG(:)  
            RDHM356_struc(n)%rdhm356(t)%NEXLAG     = tmp_NEXLAG    
            RDHM356_struc(n)%rdhm356(t)%TA_PREV    = tmp_TA_PREV   
 
           ! save output variables from local variables to global variables
            RDHM356_struc(n)%rdhm356(t)%SWINT(:)     = tmp_SWINT(:)    
            RDHM356_struc(n)%rdhm356(t)%SWHINT(:)    = tmp_SWHINT(:)   
            RDHM356_struc(n)%rdhm356(t)%TSINT(:)     = tmp_TSINT(:)    
            RDHM356_struc(n)%rdhm356(t)%FRZDUP       = tmp_FRZDUP    
            RDHM356_struc(n)%rdhm356(t)%FRZDBT       = tmp_FRZDBT    
            RDHM356_struc(n)%rdhm356(t)%FROST        = tmp_FROST    
            RDHM356_struc(n)%rdhm356(t)%ALBEDO       = tmp_ALBEDO    
            RDHM356_struc(n)%rdhm356(t)%SURF         = tmp_SURF    
            RDHM356_struc(n)%rdhm356(t)%GRND         = tmp_GRND    
            RDHM356_struc(n)%rdhm356(t)%TET          = tmp_TET    
            RDHM356_struc(n)%rdhm356(t)%EDMND        = tmp_EDMND    
            RDHM356_struc(n)%rdhm356(t)%CH           = tmp_CH    
            RDHM356_struc(n)%rdhm356(t)%CM           = tmp_CM    
            RDHM356_struc(n)%rdhm356(t)%SWE          = tmp_SWE    
            RDHM356_struc(n)%rdhm356(t)%SnowFrac     = tmp_SnowFrac    
            RDHM356_struc(n)%rdhm356(t)%SnowDepth    = tmp_SnowDepth   
            RDHM356_struc(n)%rdhm356(t)%RM           = tmp_RM    

            SoilTemp(1) = tmp_TS0 + 273.16
            SoilTemp(2) = tmp_TS1 + 273.16
            SoilTemp(3) = tmp_TS2 + 273.16
            SoilTemp(4) = tmp_TS3 + 273.16
            SoilTemp(5) = tmp_TS4 + 273.16
            SoilTemp(6) = -9999
            
            ![ 1] output variable: UZTWC (unit=mm). *** upper zone tension water storage content
            call LIS_diagnoseSurfaceOutputVar(n, t, &
                 LIS_MOC_SACUZTWC, value = RDHM356_struc(n)%rdhm356(t)%UZTWC, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, &
                 LIS_MOC_SACUZTWC, &
                 value = RDHM356_struc(n)%rdhm356(t)%UZTWC/RDHM356_struc(n)%rdhm356(t)%UZTWM, &
                 vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 2] output variable: UZFWC (unit=mm). *** upper zone free water storage content
            call LIS_diagnoseSurfaceOutputVar(n, t, &
                 LIS_MOC_SACUZFWC, value = RDHM356_struc(n)%rdhm356(t)%UZFWC, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, &
                 LIS_MOC_SACUZFWC, value = &
                 RDHM356_struc(n)%rdhm356(t)%UZFWC/RDHM356_struc(n)%rdhm356(t)%UZFWM, &
                 vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 3] output variable: LZTWC (unit=mm). *** lower zone tension water storage content
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACLZTWC, &
                 value = RDHM356_struc(n)%rdhm356(t)%LZTWC, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACLZTWC, &
                 value = RDHM356_struc(n)%rdhm356(t)%LZTWC/RDHM356_struc(n)%rdhm356(t)%LZTWM, &
                 vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 4] output variable: LZFPC (unit=mm). *** lower zone primary free water storage content
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACLZFPC, &
                 value = RDHM356_struc(n)%rdhm356(t)%LZFPC, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACLZFPC, &
                 value = RDHM356_struc(n)%rdhm356(t)%LZFPC/RDHM356_struc(n)%rdhm356(t)%LZFPM, &
                 vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 5] output variable: LZFSC (unit=mm). *** lower zone supplemental free water storage content
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACLZFSC, &
                 value = RDHM356_struc(n)%rdhm356(t)%LZFSC, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACLZFSC, &
                 value = RDHM356_struc(n)%rdhm356(t)%LZFSC/RDHM356_struc(n)%rdhm356(t)%LZFSM, &
                 vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 6] output variable: ADIMC (unit=mm). *** additional impervious area content
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACADIMPC, &
                 value = RDHM356_struc(n)%rdhm356(t)%ADIMC, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACADIMPC, &
                 value = RDHM356_struc(n)%rdhm356(t)%ADIMC &
                 / (RDHM356_struc(n)%rdhm356(t)%UZTWM +  RDHM356_struc(n)%rdhm356(t)%LZTWM), &
                 vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
!            ![ 7] output variable: TS0 (unit=K). *** first soil layer temperature (unit is Celsius in model)
!            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACHTET_TS0, value = RDHM356_struc(n)%rdhm356(t)%TS0 + 273.16, &
!                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
!           
!            ![ 8] output variable: TS1 (unit=K). *** second soil layer temperature (unit is Celsius in model)
!            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACHTET_TS1, value = RDHM356_struc(n)%rdhm356(t)%TS1 + 273.16, &
!                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
!            
!            ![ 9] output variable: TS2 (unit=K). *** third soil layer temperature (unit is Celsius in model)
!            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACHTET_TS2, value = RDHM356_struc(n)%rdhm356(t)%TS2 + 273.16, &
!                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
!            
!            ![ 10] output variable: TS3 (unit=K). *** fourth soil layer temperature (unit is Celsius in model)
!            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACHTET_TS3, value = RDHM356_struc(n)%rdhm356(t)%TS3 + 273.16, &
!                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
!            
!            ![ 11] output variable: TS4 (unit=K). *** fifth soil layer temperature (unit is Celsius in model)
!            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACHTET_TS4, value = RDHM356_struc(n)%rdhm356(t)%TS4 + 273.16, &
!                                              vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)


            ![ 12] output variable: UZTWH (unit=mm). *** unfrozen upper zone tension water
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACUZTWH, &
                 value = RDHM356_struc(n)%rdhm356(t)%UZTWH, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 13] output variable: UZFWH (unit=mm). *** unfrozen uppeer zone free water
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACUZFWH, &
                 value = RDHM356_struc(n)%rdhm356(t)%UZFWH, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 14] output variable: LZTWH (unit=mm). *** unfrozen lower zone tension water
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACLZTWH, &
                 value = RDHM356_struc(n)%rdhm356(t)%LZTWH, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 15] output variable: LZFSH (unit=mm). *** unfrozen lower zone supplemental free water
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACLZFSH, &
                 value = RDHM356_struc(n)%rdhm356(t)%LZFSH, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 16] output variable: LZFPH (unit=mm). *** unfrozen lower zone primary free water
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACLZFPH, &
                 value = RDHM356_struc(n)%rdhm356(t)%LZFPH, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 17] output variable: SMC (unit=m^3 m-3). *** volumetric content of total soil moisture at each layer 
            do i=1, 6
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST,&
                     value = RDHM356_struc(n)%rdhm356(t)%SMC(i), &
                     vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
             end do
            
            ![ 18] output variable: SH2O (unit=m^3 m-3). *** volumetric content of liquid soil moisture at each layer
            do i=1, 6
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SMLIQFRAC,&
                     value = RDHM356_struc(n)%rdhm356(t)%SH2O(i), &
                     vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
             end do
            
            ![ 19] output variable: WE (unit=kg m-2). *** snow water equivalent without liquid water
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17SWE, &
                 value = RDHM356_struc(n)%rdhm356(t)%WE, &
                 vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 20] output variable: LIQW (unit=kg m-2). *** liquid water in snow
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17LIQW, &
                 value = RDHM356_struc(n)%rdhm356(t)%LIQW, &
                 vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 21] output variable: NEGHS (unit=mm). *** negative snow heat
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17NEGHS, &
                 value = RDHM356_struc(n)%rdhm356(t)%NEGHS, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 22] output variable: TINDEX (unit=K). *** antecedent temperature index (unit is Celsius in model)
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17TINDEX, value = RDHM356_struc(n)%rdhm356(t)%TINDEX + 273.16, &
            !                                  vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 23] output variable: ACCMAX (unit=mm). *** cumulated snow water including liquid
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17ACCMAX, &
                 value = RDHM356_struc(n)%rdhm356(t)%ACCMAX, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 24] output variable: SNDPT (unit=cm). *** snow depth
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17SNDPT, value = RDHM356_struc(n)%rdhm356(t)%SNDPT, &
            !                                  vlevel=1, unit="cm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 25] output variable: SNTMP (unit=K). *** average snow temperature (unit is Celsius in model)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWT, &
                 value = RDHM356_struc(n)%rdhm356(t)%SNTMP + 273.16, &
                 vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 26] output variable: SB (unit=K). *** the last highest snow water equivalent before any snow fall (unit is Celsius in model)
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17SB, value = RDHM356_struc(n)%rdhm356(t)%SB + 273.16, &
            !                                  vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 27] output variable: SBAESC (unit=-). *** internal snow state during melt & new snow fall (checked with Victor)
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17SBAESC, value = RDHM356_struc(n)%rdhm356(t)%SBAESC, &
            !                                  vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 28] output variable: SBWS (unit=-). *** internal snow state during melt & new snow fall (checked with Victor)
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17SBWS, value = RDHM356_struc(n)%rdhm356(t)%SBWS, &
            !                                  vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 29] output variable: STORAGE (unit=mm). *** snow liquid water attenuation storage
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17STORAGE, value = RDHM356_struc(n)%rdhm356(t)%STORAGE, &
            !                                  vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 30] output variable: AEADJ (unit=-). *** adjusted areal snow cover fraction
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17AEADJ, &
                 value = RDHM356_struc(n)%rdhm356(t)%AEADJ, &
                 vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 31] output variable: SWINT (unit=-). *** total volumetric soil moisture contents at desired soil layers (can be different from soil layers)
            do i=1, RDHM356_struc(n)%NDINTW
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACSWINT, &
                     value = RDHM356_struc(n)%rdhm356(t)%SWINT(i), &
                     vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 32] output variable: SWHINT (unit=-). *** liquid volumetric soil moisture contents at desired soil layers (can be different from soil layers)
            do i=1, RDHM356_struc(n)%NDINTW
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACSWHINT, &
                     value = RDHM356_struc(n)%rdhm356(t)%SWHINT(i), &
                     vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 33] output variable: TSINT (unit=-). *** soil temperature at desired soil layers (can be different from soil layers)
            do i=1, RDHM356_struc(n)%NDSINT
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACTSINT,&
                     value = RDHM356_struc(n)%rdhm356(t)%TSINT(i), &
                     vlevel=i, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 34] output variable: FRZDUP (unit=m). *** depth of the upper border of frozen ground from surface  
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACHTET_FRZDUP, value = RDHM356_struc(n)%rdhm356(t)%FRZDUP, &
            !                                  vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 35] output variable: FRZDBT (unit=m). *** depth of the bottom border of frozen ground from surface 
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACHTET_FRZDBT, value = RDHM356_struc(n)%rdhm356(t)%FRZDBT, &
            !                                  vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 36] output variable: FROST (unit=-). *** frost index
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SACFROST, &
                 value = RDHM356_struc(n)%rdhm356(t)%FROST, &
                 vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 37] output variable: ALBEDO (unit=-). *** land surface albedo 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ALBEDO, &
                 value = RDHM356_struc(n)%rdhm356(t)%ALBEDO, &
                 vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 38] output variable: SURF (unit=kg m-2 s-1). *** Qs <=> SURF simulated fast runoff (surface runoff)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QS, &
                 value = RDHM356_struc(n)%rdhm356(t)%SURF, &
                 vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
            
            ![ 39] output variable: GRND (unit=kg m-2 s-1). *** Qsb <=> GRND simulated slow runoff (baseflow)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QSB, &
                 value = RDHM356_struc(n)%rdhm356(t)%GRND, &
                 vlevel=1, unit="kg m-2 s-1", direction="OUT", surface_type = LIS_rc%lsm_index)
            
            ![ 40] output variable: TET (unit=kg m-2 s-1). *** Evap <=> TET simulated actual evapotranspiration
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVAP, &
                 value = 0.01 * RDHM356_struc(n)%rdhm356(t)%TET, &
                 vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_EVAP, &
                 value = 3600 * RDHM356_struc(n)%rdhm356(t)%TET, &
                 vlevel=1, unit="mm hr-1", direction="UP", surface_type = LIS_rc%lsm_index)
            
            ![ 41] output variable: EDMND (unit=kg m-2 s-1). *** PotEvap <=> EDMND potential evapotranspiration  
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_POTEVAP,&
                 value = 0.01 * RDHM356_struc(n)%rdhm356(t)%EDMND, &
                 vlevel=1, unit="kg m-2 s-1", direction="UP", surface_type = LIS_rc%lsm_index)
            
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_POTEVAP, &
                 value = 3600 * RDHM356_struc(n)%rdhm356(t)%EDMND, &
                 vlevel=1, unit="mm hr-1", direction="UP", surface_type = LIS_rc%lsm_index)
 
            ![ 42] output variable: CH (unit=s m-1). *** surface layer exchage coefficient for heat and moisture
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CH, &
                 value = RDHM356_struc(n)%rdhm356(t)%CH, &
                 vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 43] output variable: CM (unit=s m-1). *** surface layer exchange coefficient for momentum (drag coefficient)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CM, &
                 value = RDHM356_struc(n)%rdhm356(t)%CM, &
                 vlevel=1, unit="m s-1", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 44] output variable: SnowFrac (unit=-). *** snow cover fraction, Snow-17
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWCOVER, &
                 value = RDHM356_struc(n)%rdhm356(t)%SnowFrac, &
                 vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 45] output variable: SWE (unit=kg m-2). *** snow water equivalent, Snow-17
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SWE, &
                 value = RDHM356_struc(n)%rdhm356(t)%SWE, &
                 vlevel=1, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 46] output variable: SnowDepth (unit=m). *** snow depth, Snow-17
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDEPTH, &
                 value = RDHM356_struc(n)%rdhm356(t)%SnowDepth, &
                 vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDEPTH, &
                 value = RDHM356_struc(n)%rdhm356(t)%SnowDepth * 100.0, &
                 vlevel=1, unit="cm", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOWDEPTH, &
                 value = RDHM356_struc(n)%rdhm356(t)%SnowDepth * 1000.0, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 47] output variable: RM (unit=mm). *** rain + melt
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SNOW17RMLT, &
                 value = RDHM356_struc(n)%rdhm356(t)%RM, &
                 vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![48] soil temperature K
            do i=1, 6
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILTEMP, &
                     value = SoilTemp(i),                       &
                     vlevel=i, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ! reset forcing variables to zeros
            RDHM356_struc(n)%rdhm356(t)%Tair = 0.0
            RDHM356_struc(n)%rdhm356(t)%Psurf = 0.0
            RDHM356_struc(n)%rdhm356(t)%Wind_E = 0.0
            RDHM356_struc(n)%rdhm356(t)%Wind_N = 0.0
            RDHM356_struc(n)%rdhm356(t)%Qair = 0.0
            RDHM356_struc(n)%rdhm356(t)%Rainf = 0.0
            RDHM356_struc(n)%rdhm356(t)%Snowf = 0.0
            RDHM356_struc(n)%rdhm356(t)%Swdown = 0.0
            RDHM356_struc(n)%rdhm356(t)%Lwdown = 0.0
        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        RDHM356_struc(n)%forc_count = 0 
    endif ! end of alarmCheck loop 
    
    deallocate( tmp_PET_MON )
    deallocate( tmp_PETADJ_MON )
    deallocate( tmp_GRN_MON )
    deallocate( tmp_SMC )
    deallocate( tmp_SH2O )
    deallocate( tmp_DSINTW )
    deallocate( tmp_DSINT )
    deallocate( tmp_SWINT )
    deallocate( tmp_SWHINT )
    deallocate( tmp_TSINT )
    deallocate( tmp_ADC )
    deallocate( tmp_EXLAG )
end subroutine RDHM356_main
