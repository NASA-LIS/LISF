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
! !ROUTINE: AWRAL600_main
! \label{AWRAL600_main}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   12/18/18: Wendy Sharples, Shugong Wang; initial implementation for AWRAL600 with LIS-7
!
! !INTERFACE:
subroutine AWRAL600_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod 
    use AWRAL600_lsmMod
   !use other modules
  
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
!  This is the entry point for calling the AWRAL600 physics.
!  This routine calls the {\tt awral_driver_600 } routine that performs the
!  land surface computations, to solve for water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for AWRAL600
    integer, parameter            :: dp = selected_real_kind(15, 307) ! to store tmp vars as doubles - convert back
    integer                       :: tmp_n                  ! nest id [-]
    real                          :: tmp_latitude           ! latitude in decimal degree [rad]
    real                          :: tmp_longitude           ! longitude in decimal year [rad]
    integer                       :: tmp_year               ! year of the currrent time step [-]
    integer                       :: tmp_month              ! month of the current time step [-]
    integer                       :: tmp_day                ! day of the current time step [-]
    integer                       :: tmp_hour               ! hour of the current time step [-]
    integer                       :: tmp_minute             ! minute of the current time step [-]
    real(kind=dp)                 :: tmp_Tair               ! average air temperature [K]
    real(kind=dp)                 :: tmp_Swdown             ! downward shortwave radiation [W/m2]
    real(kind=dp)                 :: tmp_Rainf              ! daily gross precipitation [kg/m2s]
    real(kind=dp)                 :: tmp_Qair               ! actual vapour pressure [kg/kg]
    real(kind=dp)                 :: tmp_Wind_E             ! 2m wind magnitude [m/s]
    real(kind=dp)                 :: tmp_Swdirect           ! expected downwelling shortwave radiation on a cloudless day [W/m2]
    real(kind=dp)                 :: tmp_e0                 ! potential evaporation [mm/d]
    real(kind=dp)                 :: tmp_etot               ! actual evapotranspiration [mm/d]
    real(kind=dp)                 :: tmp_dd                 ! vertical drainage from the bottom of the deep soil layer [mm]
    real(kind=dp)                 :: tmp_s0_avg             ! water storage in the surface soil layer [mm]
    real(kind=dp)                 :: tmp_ss_avg             ! water content of the shallow soil store [mm]
    real(kind=dp)                 :: tmp_sd_avg             ! water content of the deep soil store [mm]
    real(kind=dp)                 :: tmp_qtot               ! total discharge to stream [mm]
    real(kind=dp)                 :: tmp_sr                 ! volume of water in the surface water store [mm]
    real(kind=dp)                 :: tmp_sg                 ! groundwater storage in the unconfined aquifer [mm]
    real(kind=dp), allocatable    :: tmp_s0(:)              ! water storage in the surface soil layer for each hru [mm]
    real(kind=dp), allocatable    :: tmp_ss(:)              ! water content of the shallow soil store for each hru [mm]
    real(kind=dp), allocatable    :: tmp_sd(:)              ! water content of the deep soil store for each hru [mm]
    real(kind=dp), allocatable    :: tmp_mleaf(:)           ! leaf biomass [kg/m2]
    integer                       :: tmp_nhru               ! number of hydrologic response units [-]
    integer                       :: tmp_nhypsbins          ! number of hypsometric curve percentile bins [-]
    real(kind=dp)                 :: tmp_slope_coeff        ! scaling factor for slope [-]
    real(kind=dp)                 :: tmp_pair               ! air pressure [Pa]
    real(kind=dp)                 :: tmp_kr_coeff           ! scaling factor for ratio of saturated hydraulic conductivity [-]
    real(kind=dp)                 :: tmp_k_rout             ! rate coefficient controlling discharge to stream [-]
    real(kind=dp)                 :: tmp_kssat              ! saturated hydraulic conductivity of shallow soil layer [mm/d]
    real(kind=dp)                 :: tmp_prefr              ! reference value for precipitation [mm]
    real(kind=dp)                 :: tmp_s0max              ! maximum storage of the surface soil layer [mm]
    real(kind=dp)                 :: tmp_slope              ! slope of the land surface [%]
    real(kind=dp)                 :: tmp_ssmax              ! maximum storage of the shallow soil layer [mm]
    real(kind=dp)                 :: tmp_k_gw               ! groundwater drainage coefficient [1/d]
    real(kind=dp)                 :: tmp_kr_sd              ! routing delay factor for the deep layer [-]
    real(kind=dp)                 :: tmp_kr_0s              ! routing delay factor for the surface layer [-]
    real(kind=dp)                 :: tmp_k0sat              ! saturated hydraulic conductivity of surface soil layer [mm/d]
    real(kind=dp)                 :: tmp_sdmax              ! maximum storage of the deep soil layer [mm]
    real(kind=dp)                 :: tmp_kdsat              ! saturated hydraulic conductivity of shallow soil layer [mm/d]
    real(kind=dp)                 :: tmp_ne                 ! effective porosity [-]
    real(kind=dp), allocatable    :: tmp_height(:)             ! elevation of a point on the hypsometric curve [m]
    real(kind=dp), allocatable    :: tmp_hypsperc(:)           ! hypsometric curve distribution percentile bins [%]
    real(kind=dp), allocatable    :: tmp_alb_dry(:)         ! dry soil albedo for each hru [-]
    real(kind=dp), allocatable    :: tmp_alb_wet(:)         ! wet soil albedo for each hru [-]
    real(kind=dp), allocatable    :: tmp_cgsmax(:)          ! coefficient relating vegetation photosynthetic capacity to maximum stomatal conductance for each hru [m/s]
    real(kind=dp), allocatable    :: tmp_er_frac_ref(:)     ! specific ratio of the mean evaporation rate and the mean rainfall intensity during storms for each hru [-]
    real(kind=dp), allocatable    :: tmp_fsoilemax(:)       ! soil evaporation scaling factor corresponding to unlimited soil water supply for each hru [-]
    real(kind=dp), allocatable    :: tmp_lairef(:)          ! reference leaf area index (at which fv = 0.63) for each hru [-]
    real(kind=dp), allocatable    :: tmp_rd(:)              ! rooting depth for each hru [m]
    real(kind=dp), allocatable    :: tmp_s_sls(:)           ! specific canopy rainfall storage per unit leaf area for each hru [mm]
    real(kind=dp), allocatable    :: tmp_sla(:)             ! specific leaf area for each hru [m2/kg]
    real(kind=dp), allocatable    :: tmp_tgrow(:)           ! characteristic time scale for vegetation growth towards equilibrium for each hru [d]
    real(kind=dp), allocatable    :: tmp_tsenc(:)           ! characteristic time scale for vegetation senescence towards equilibrium for each hru [d]
    real(kind=dp), allocatable    :: tmp_ud0(:)             ! maximum possible root water uptake from the deep soil store for each hru [mm/d]
    real(kind=dp), allocatable    :: tmp_us0(:)             ! maximum possible root water uptake from the shallow soil store for each hru [mm/d]
    real(kind=dp), allocatable    :: tmp_vc(:)              ! vegetation photosynthetic capacity index per unit canopy cover for each hru [-]
    real(kind=dp), allocatable    :: tmp_w0lime(:)          ! limiting the value of the relative soil moisture content of the top soil layer at which evaporation is reduced for each hru [-]
    real(kind=dp), allocatable    :: tmp_w0ref_alb(:)       ! Reference value of w0 that determines the rate of albedo decrease with wetness for each hru [-]
    real(kind=dp), allocatable    :: tmp_wdlimu(:)          ! water-limiting relative water content of the deep soil store for each hru [-]
    real(kind=dp), allocatable    :: tmp_wslimu(:)          ! water-limiting relative water content of the shallow soil store for each hru [-]
    real(kind=dp), allocatable    :: tmp_fhru(:)            ! fraction of the cell which contains shallow and deep rooted vegetation [-]
    real(kind=dp), allocatable    :: tmp_hveg(:)            ! vegetation height for each hru [-]
    real(kind=dp), allocatable    :: tmp_laimax(:)          ! leaf area index max for each hru [-]
    integer                       :: tmp_timesteps          ! number of daily timesteps [-]

    allocate( tmp_hypsperc( AWRAL600_struc(n)%nhypsbins ) )
    allocate( tmp_height( AWRAL600_struc(n)%nhypsbins ) )
    allocate( tmp_s0( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_ss( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_sd( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_mleaf( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_alb_dry( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_alb_wet( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_cgsmax( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_er_frac_ref( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_fsoilemax( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_lairef( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_rd( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_s_sls( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_sla( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_tgrow( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_tsenc( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_ud0( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_us0( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_vc( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_w0lime( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_w0ref_alb( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_wdlimu( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_wslimu( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_fhru( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_hveg( AWRAL600_struc(n)%nhru ) )
    allocate( tmp_laimax( AWRAL600_struc(n)%nhru ) )

    ! check AWRAL600 alarm. If alarm is ring, run model. 
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "AWRAL600 model alarm")
    if (alarmCheck) Then
        do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            dt = LIS_rc%ts
            row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
            lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

            ! retrieve forcing data from AWRAL600_struc(n)%awral600(t) and assign to local variables
            ! Tair: average air temperature
            tmp_Tair         = AWRAL600_struc(n)%awral600(t)%Tair    / AWRAL600_struc(n)%forc_count
            tmp_Tair = ANINT(1000000.0*tmp_Tair)/1000000.0
 
            ! Swdown: downward shortwave radiation
            tmp_Swdown       = 0.0 + AWRAL600_struc(n)%awral600(t)%Swdown   / AWRAL600_struc(n)%forc_count
            tmp_Swdown = ANINT(1000000.0*tmp_Swdown)/1000000.0
 
            ! Rainf: daily gross precipitation
            tmp_Rainf        = 0.0 + AWRAL600_struc(n)%awral600(t)%Rainf    / AWRAL600_struc(n)%forc_count
	    tmp_Rainf = ANINT(1000000.0*tmp_Rainf)/1000000.0 

            ! Qair: actual vapour pressure
            tmp_Qair         = 0.0 + AWRAL600_struc(n)%awral600(t)%Qair     / AWRAL600_struc(n)%forc_count
            tmp_Qair = ANINT(1000000.0*tmp_Qair)/1000000.0 

            ! Wind_E: 2m wind magnitude
            tmp_Wind_E       = 0.0 + AWRAL600_struc(n)%awral600(t)%Wind_E   / AWRAL600_struc(n)%forc_count
	    tmp_Wind_E = ANINT(1000000.0*tmp_Wind_E)/1000000.0 

            ! Swdirect: expected downwelling shortwave radiation on a cloudless day
            tmp_Swdirect     = 0.0 + AWRAL600_struc(n)%awral600(t)%Swdirect / AWRAL600_struc(n)%forc_count
	    tmp_Swdirect = ANINT(1000000.0*tmp_Swdirect)/1000000.0

            ! check validity of Tair
            if(tmp_Tair .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable Tair in AWRAL600"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Swdown
            if(tmp_Swdown .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable Swdown in AWRAL600"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Rainf
            if(tmp_Rainf .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable Rainf in AWRAL600"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Qair
            if(tmp_Qair .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable Qair in AWRAL600"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Wind_E
            if(tmp_Wind_E .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable Wind_E in AWRAL600"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            ! check validity of Swdirect
            if(tmp_Swdirect .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "[ERR] undefined value found for forcing variable Swdirect in AWRAL600"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif
            
            tmp_latitude  = lat
            tmp_longitude = lon
            tmp_year   = LIS_rc%yr
            tmp_month  = LIS_rc%mo
            tmp_day    = LIS_rc%da
            tmp_hour   = LIS_rc%hr
            tmp_minute = LIS_rc%mn
 
            ! get parameters - single precision 
            tmp_nhru                                = AWRAL600_struc(n)%nhru
            tmp_nhypsbins                           = AWRAL600_struc(n)%nhypsbins
            tmp_slope_coeff                         = AWRAL600_struc(n)%slope_coeff 
	    tmp_slope_coeff                	    = ANINT(1000000.0*tmp_slope_coeff)/1000000.0
            tmp_pair                                = AWRAL600_struc(n)%pair 
	    tmp_pair				    = ANINT(1000000.0*tmp_pair)/1000000.0                           
            tmp_kr_coeff                            = AWRAL600_struc(n)%kr_coeff
	    tmp_kr_coeff                            = ANINT(1000000.0*tmp_kr_coeff)/1000000.0                        
            tmp_k_rout                              = AWRAL600_struc(n)%awral600(t)%k_rout
	    tmp_k_rout                              = ANINT(1000000.0*tmp_k_rout)/1000000.0                          
            tmp_kssat                               = AWRAL600_struc(n)%awral600(t)%kssat                         
	    tmp_kssat                               = ANINT(1000000.0*tmp_kssat)/1000000.0
            tmp_prefr                               = AWRAL600_struc(n)%awral600(t)%prefr                          
	    tmp_prefr                               = ANINT(1000000.0*tmp_prefr)/1000000.0
            tmp_s0max                               = AWRAL600_struc(n)%awral600(t)%s0max                           
	    tmp_s0max                               = ANINT(1000000.0*tmp_s0max)/1000000.0
            tmp_slope                               = AWRAL600_struc(n)%awral600(t)%slope                           
	    tmp_slope                               = ANINT(1000000.0*tmp_slope)/1000000.0
            tmp_ssmax                               = AWRAL600_struc(n)%awral600(t)%ssmax                          
	    tmp_ssmax                               = ANINT(1000000.0*tmp_ssmax)/1000000.0
            tmp_k_gw                                = AWRAL600_struc(n)%awral600(t)%k_gw                            
	    tmp_k_gw                                = ANINT(1000000.0*tmp_k_gw)/1000000.0
            tmp_kr_sd                               = AWRAL600_struc(n)%awral600(t)%kr_sd                           
	    tmp_kr_sd                               = ANINT(1000000.0*tmp_kr_sd)/1000000.0
            tmp_kr_0s                               = AWRAL600_struc(n)%awral600(t)%kr_0s                           
	    tmp_kr_0s                               = ANINT(1000000.0*tmp_kr_0s)/1000000.0
            tmp_k0sat                               = AWRAL600_struc(n)%awral600(t)%k0sat                           
	    tmp_k0sat                               = ANINT(1000000.0*tmp_k0sat)/1000000.0
            tmp_sdmax                               = AWRAL600_struc(n)%awral600(t)%sdmax                           
	    tmp_sdmax                               = ANINT(1000000.0*tmp_sdmax)/1000000.0
            tmp_kdsat                               = AWRAL600_struc(n)%awral600(t)%kdsat                           
	    tmp_kdsat                               = ANINT(1000000.0*tmp_kdsat)/1000000.0
            tmp_ne                                  = AWRAL600_struc(n)%awral600(t)%ne                              
	    tmp_ne                                  = ANINT(1000000.0*tmp_ne)/1000000.0
            tmp_height(:)                           = AWRAL600_struc(n)%awral600(t)%height(:)                        
	    tmp_height(:)                           = ANINT(1000000.0*tmp_height(:))/1000000.0
            tmp_hypsperc(:)                         = AWRAL600_struc(n)%hypsperc(:)                        
	    tmp_hypsperc(:)                         = ANINT(1000000.0*tmp_hypsperc(:))/1000000.0
            tmp_alb_dry(:)                          = AWRAL600_struc(n)%alb_dry(:)                        
	    tmp_alb_dry(:)                          = ANINT(1000000.0*tmp_alb_dry(:))/1000000.0
            tmp_alb_wet(:)                          = AWRAL600_struc(n)%alb_wet(:)                         
	    tmp_alb_wet(:)                          = ANINT(1000000.0*tmp_alb_wet(:))/1000000.0
            tmp_cgsmax(:)                           = AWRAL600_struc(n)%cgsmax(:)                          
	    tmp_cgsmax(:)                           = ANINT(1000000.0*tmp_cgsmax(:))/1000000.0
            tmp_er_frac_ref(:)                      = AWRAL600_struc(n)%er_frac_ref(:)                    
	    tmp_er_frac_ref(:)                      = ANINT(1000000.0*tmp_er_frac_ref(:))/1000000.0
            tmp_fsoilemax(:)                        = AWRAL600_struc(n)%fsoilemax(:)                       
	    tmp_fsoilemax(:)                        = ANINT(1000000.0*tmp_fsoilemax(:))/1000000.0
            tmp_lairef(:)                           = AWRAL600_struc(n)%lairef(:)                          
	    tmp_lairef(:)                           = ANINT(1000000.0*tmp_lairef(:))/1000000.0
            tmp_rd(:)                               = AWRAL600_struc(n)%rd(:)                              
	    tmp_rd(:)                               = ANINT(1000000.0*tmp_rd(:))/1000000.0
            tmp_s_sls(:)                            = AWRAL600_struc(n)%s_sls(:)                           
	    tmp_s_sls(:)                            = ANINT(1000000.0*tmp_s_sls(:))/1000000.0
            tmp_sla(:)                              = AWRAL600_struc(n)%sla(:)                             
	    tmp_sla(:)                              = ANINT(1000000.0*tmp_sla(:))/1000000.0
            tmp_tgrow(:)                            = AWRAL600_struc(n)%tgrow(:)                           
	    tmp_tgrow(:)                            = ANINT(1000000.0*tmp_tgrow(:))/1000000.0
            tmp_tsenc(:)                            = AWRAL600_struc(n)%tsenc(:)                           
	    tmp_tsenc(:)                            = ANINT(1000000.0*tmp_tsenc(:))/1000000.0
            tmp_ud0(:)                              = AWRAL600_struc(n)%ud0(:)                            
	    tmp_ud0(:)                              = ANINT(1000000.0*tmp_ud0(:))/1000000.0
            tmp_us0(:)                              = AWRAL600_struc(n)%us0(:)                             
	    tmp_us0(:)                              = ANINT(1000000.0*tmp_us0(:))/1000000.0
            tmp_vc(:)                               = AWRAL600_struc(n)%vc(:)                              
	    tmp_vc(:)                               = ANINT(1000000.0*tmp_vc(:))/1000000.0
            tmp_w0lime(:)                           = AWRAL600_struc(n)%w0lime(:)                          
	    tmp_w0lime(:)                           = ANINT(1000000.0*tmp_w0lime(:))/1000000.0
            tmp_w0ref_alb(:)                        = AWRAL600_struc(n)%w0ref_alb(:)                      
	    tmp_w0ref_alb(:)                        = ANINT(1000000.0*tmp_w0ref_alb(:))/1000000.0
            tmp_wdlimu(:)                           = AWRAL600_struc(n)%wdlimu(:)                          
	    tmp_wdlimu(:)                           = ANINT(1000000.0*tmp_wdlimu(:))/1000000.0
            tmp_wslimu(:)                           = AWRAL600_struc(n)%wslimu(:)                         
	    tmp_wslimu(:)                           = ANINT(1000000.0*tmp_wslimu(:))/1000000.0
            tmp_fhru(:)                             = AWRAL600_struc(n)%awral600(t)%fhru(:)                            
	    tmp_fhru(:)                             = ANINT(1000000.0*tmp_fhru(:))/1000000.0
            tmp_hveg(:)                             = AWRAL600_struc(n)%awral600(t)%hveg(:)                            
	    tmp_hveg(:)                             = ANINT(1000000.0*tmp_hveg(:))/1000000.0
            tmp_laimax(:)                           = AWRAL600_struc(n)%awral600(t)%laimax(:)                          
	    tmp_laimax(:)                           = ANINT(1000000.0*tmp_laimax(:))/1000000.0
            tmp_timesteps                           = AWRAL600_struc(n)%timesteps
 
            ! get state variables
            tmp_sr       = ANINT(1000000.0*AWRAL600_struc(n)%awral600(t)%sr)/1000000.0   
            tmp_sg       = ANINT(1000000.0*AWRAL600_struc(n)%awral600(t)%sg)/1000000.0   
            tmp_s0(:)    = ANINT(1000000.0*AWRAL600_struc(n)%awral600(t)%s0(:))/1000000.0   
            tmp_ss(:)    = ANINT(1000000.0*AWRAL600_struc(n)%awral600(t)%ss(:))/1000000.0   
            tmp_sd(:)    = ANINT(1000000.0*AWRAL600_struc(n)%awral600(t)%sd(:))/1000000.0   
            tmp_mleaf(:) = ANINT(1000000.0*AWRAL600_struc(n)%awral600(t)%mleaf(:))/1000000.0

            ! call model physics 
            call awral_driver_600(tmp_Tair              , & ! in    - average air temperature [K]
                                  tmp_Swdown            , & ! in    - downward shortwave radiation [W/m2]
                                  tmp_Rainf             , & ! in    - daily gross precipitation [kg/m2s]
                                  tmp_Qair              , & ! in    - actual vapour pressure [kg/kg]
                                  tmp_Wind_E            , & ! in    - 2m wind magnitude [m/s]
                                  tmp_Swdirect          , & ! in    - expected downwelling shortwave radiation on a cloudless day [W/m2]
                                  tmp_e0                , & ! out   - potential evaporation [mm/d]
                                  tmp_etot              , & ! out   - actual evapotranspiration [mm/d]
                                  tmp_dd                , & ! out   - vertical drainage from the bottom of the deep soil layer [mm]
                                  tmp_s0_avg            , & ! out   - water storage in the surface soil layer [mm]
                                  tmp_ss_avg            , & ! out   - water content of the shallow soil store [mm]
                                  tmp_sd_avg            , & ! out   - water content of the deep soil store [mm]
                                  tmp_qtot              , & ! out   - total discharge to stream [mm]
                                  tmp_sr                , & ! inout - volume of water in the surface water store [mm]
                                  tmp_sg                , & ! inout - groundwater storage in the unconfined aquifer [mm]
                                  tmp_s0                , & ! inout - water storage in the surface soil layer for each hru [mm]
                                  tmp_ss                , & ! inout - water content of the shallow soil store for each hru [mm]
                                  tmp_sd                , & ! inout - water content of the deep soil store for each hru [mm]
                                  tmp_mleaf             , & ! inout - leaf biomass [kg/m2]
                                  !not included but could be tmp_nhru              , & ! in    - number of hydrologic response units [-]
                                  !not included but could be tmp_nhypsbins         , & ! in    - number of hypsometric percentile bins [-]
                                  tmp_slope_coeff       , & ! in    - scaling factor for slope [-]
                                  tmp_pair              , & ! in    - air pressure [Pa]
                                  tmp_kr_coeff          , & ! in    - scaling factor for ratio of saturated hydraulic conductivity [-]
                                  tmp_k_rout            , & ! in    - rate coefficient controlling discharge to stream [-]
                                  tmp_kssat             , & ! in    - saturated hydraulic conductivity of shallow soil layer [mm/d]
                                  tmp_prefr             , & ! in    - reference value for precipitation [mm]
                                  tmp_s0max             , & ! in    - maximum storage of the surface soil layer [mm]
                                  tmp_slope             , & ! in    - slope of the land surface [%]
                                  tmp_ssmax             , & ! in    - maximum storage of the shallow soil layer [mm]
                                  tmp_k_gw              , & ! in    - groundwater drainage coefficient [1/d]
                                  tmp_kr_sd             , & ! in    - routing delay factor for the deep layer [-]
                                  tmp_kr_0s             , & ! in    - routing delay factor for the surface layer [-]
                                  tmp_k0sat             , & ! in    - saturated hydraulic conductivity of surface soil layer [mm/d]
                                  tmp_sdmax             , & ! in    - maximum storage of the deep soil layer [mm]
                                  tmp_kdsat             , & ! in    - saturated hydraulic conductivity of shallow soil layer [mm/d]
                                  tmp_ne                , & ! in    - effective porosity [-]
                                  tmp_height            , & ! in    - elevation of a point on the hypsometric curve [m]
                                  tmp_hypsperc          , & ! in    - hypsometric curve distribution percentile bins [%]
                                  tmp_alb_dry           , & ! in    - dry soil albedo for each hru [-]
                                  tmp_alb_wet           , & ! in    - wet soil albedo for each hru [-]
                                  tmp_cgsmax            , & ! in    - coefficient relating vegetation photosynthetic capacity to maximum stomatal conductance for each hru [m/s]
                                  tmp_er_frac_ref       , & ! in    - specific ratio of the mean evaporation rate and the mean rainfall intensity during storms for each hru [-]
                                  tmp_fsoilemax         , & ! in    - soil evaporation scaling factor corresponding to unlimited soil water supply for each hru [-]
                                  tmp_lairef            , & ! in    - reference leaf area index (at which fv = 0.63) for each hru [-]
                                  tmp_rd                , & ! in    - rooting depth for each hru [m]
                                  tmp_s_sls             , & ! in    - specific canopy rainfall storage per unit leaf area for each hru [mm]
                                  tmp_sla               , & ! in    - specific leaf area for each hru [m2/kg]
                                  tmp_tgrow             , & ! in    - characteristic time scale for vegetation growth towards equilibrium for each hru [d]
                                  tmp_tsenc             , & ! in    - characteristic time scale for vegetation senescence towards equilibrium for each hru [d]
                                  tmp_ud0               , & ! in    - maximum possible root water uptake from the deep soil store for each hru [mm/d]
                                  tmp_us0               , & ! in    - maximum possible root water uptake from the shallow soil store for each hru [mm/d]
                                  tmp_vc                , & ! in    - vegetation photosynthetic capacity index per unit canopy cover for each hru [-]
                                  tmp_w0lime            , & ! in    - limiting the value of the relative soil moisture content of the top soil layer at which evaporation is reduced for each hru [-]
                                  tmp_w0ref_alb         , & ! in    - Reference value of w0 that determines the rate of albedo decrease with wetness for each hru [-]
                                  tmp_wdlimu            , & ! in    - water-limiting relative water content of the deep soil store for each hru [-]
                                  tmp_wslimu            , & ! in    - water-limiting relative water content of the shallow soil store for each hru [-]
                                  tmp_fhru              , & ! in    - fraction of the cell which contains shallow and deep rooted vegetation [-]
                                  tmp_hveg              , & ! in    - vegetation height for each hru [-]
                                  tmp_laimax            , & ! in    - leaf area index max for each hru [-]
                                  tmp_timesteps         ) ! in    - number of daily timesteps [-] SHOULD ALWAYS BE 1 LIS takes care of the time evolution
 
    
            ! save state variables from local variables to global variables
            AWRAL600_struc(n)%awral600(t)%sr       = SNGL(tmp_sr)      
            AWRAL600_struc(n)%awral600(t)%sg       = SNGL(tmp_sg)      
            AWRAL600_struc(n)%awral600(t)%s0(:)    = SNGL(tmp_s0(:))   
            AWRAL600_struc(n)%awral600(t)%ss(:)    = SNGL(tmp_ss(:))   
            AWRAL600_struc(n)%awral600(t)%sd(:)    = SNGL(tmp_sd(:))   
            AWRAL600_struc(n)%awral600(t)%mleaf(:) = SNGL(tmp_mleaf(:))

            ! save output variables from local variables to global variables
            AWRAL600_struc(n)%awral600(t)%e0        = SNGL(tmp_e0)       
            AWRAL600_struc(n)%awral600(t)%etot      = SNGL(tmp_etot)     
            AWRAL600_struc(n)%awral600(t)%dd        = SNGL(tmp_dd)       
            AWRAL600_struc(n)%awral600(t)%s0_avg    = SNGL(tmp_s0_avg)  
            AWRAL600_struc(n)%awral600(t)%ss_avg    = SNGL(tmp_ss_avg)   
            AWRAL600_struc(n)%awral600(t)%sd_avg    = SNGL(tmp_sd_avg)   
            AWRAL600_struc(n)%awral600(t)%qtot      = SNGL(tmp_qtot)     
            
            ![ 1] output variable: sr (unit=mm). *** volume of water in the surface water store
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SR, value = AWRAL600_struc(n)%awral600(t)%sr, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 2] output variable: sg (unit=mm). *** groundwater storage in the unconfined aquifer
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SG, value = AWRAL600_struc(n)%awral600(t)%sg, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 3] output variable: s0 (unit=mm). ***  water storage in the surface soil layer for each hru
            do i=1, AWRAL600_struc(n)%nhru
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_S0_HRU, value = AWRAL600_struc(n)%awral600(t)%s0(i), &
                                                  vlevel=i, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 4] output variable: ss (unit=mm). ***  water content of the shallow soil store for each hru
            do i=1, AWRAL600_struc(n)%nhru
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SS_HRU, value = AWRAL600_struc(n)%awral600(t)%ss(i), &
                                                  vlevel=i, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 5] output variable: sd (unit=mm). ***   water content of the deep soil store for each hru
            do i=1, AWRAL600_struc(n)%nhru
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SD_HRU, value = AWRAL600_struc(n)%awral600(t)%sd(i), &
                                                  vlevel=i, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 6] output variable: mleaf (unit=kg/m2). ***  leaf biomass for each hru
            do i=1, AWRAL600_struc(n)%nhru
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_MLEAF_HRU, value = AWRAL600_struc(n)%awral600(t)%mleaf(i), &
                                                  vlevel=i, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            
            ![ 7] output variable: e0 (unit=mm). *** potential evaporation
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_E0, value = AWRAL600_struc(n)%awral600(t)%e0, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 8] output variable: etot (unit=mm). *** actual evapotranspiration
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_ETOT, value = AWRAL600_struc(n)%awral600(t)%etot, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 9] output variable: dd (unit=mm). *** vertical drainage from the bottom of the deep soil layer
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_DD, value = AWRAL600_struc(n)%awral600(t)%dd, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 10] output variable: s0_avg (unit=mm). *** water storage in the surface soil layer
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_S0, value = AWRAL600_struc(n)%awral600(t)%s0_avg, &
                                              vlevel=1, unit="mm", direction="- ", surface_type = LIS_rc%lsm_index)
            
            ![ 11] output variable: ss_avg (unit=mm). *** water content of the shallow soil store
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SS, value = AWRAL600_struc(n)%awral600(t)%ss_avg, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 12] output variable: sd_avg (unit=mm). *** water content of the deep soil store
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SD, value = AWRAL600_struc(n)%awral600(t)%sd_avg, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            
            ![ 13] output variable: qtot (unit=mm). *** total discharge to stream
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_QTOT, value = AWRAL600_struc(n)%awral600(t)%qtot, &
                                              vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            ! reset forcing variables to zeros
            AWRAL600_struc(n)%awral600(t)%Tair = 0.0
            AWRAL600_struc(n)%awral600(t)%Swdown = 0.0
            AWRAL600_struc(n)%awral600(t)%Rainf = 0.0
            AWRAL600_struc(n)%awral600(t)%Qair = 0.0
            AWRAL600_struc(n)%awral600(t)%Wind_E = 0.0
            AWRAL600_struc(n)%awral600(t)%Swdirect = 0.0
        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        AWRAL600_struc(n)%forc_count = 0 
    endif ! end of alarmCheck loop 
    
    deallocate( tmp_s0 )
    deallocate( tmp_ss )
    deallocate( tmp_sd )
    deallocate( tmp_mleaf )
    deallocate( tmp_alb_dry )
    deallocate( tmp_alb_wet )
    deallocate( tmp_cgsmax )
    deallocate( tmp_er_frac_ref )
    deallocate( tmp_fsoilemax )
    deallocate( tmp_lairef )
    deallocate( tmp_rd )
    deallocate( tmp_s_sls )
    deallocate( tmp_sla )
    deallocate( tmp_tgrow )
    deallocate( tmp_tsenc )
    deallocate( tmp_ud0 )
    deallocate( tmp_us0 )
    deallocate( tmp_vc )
    deallocate( tmp_w0lime )
    deallocate( tmp_w0ref_alb )
    deallocate( tmp_wdlimu )
    deallocate( tmp_wslimu )
    deallocate( tmp_height )
    deallocate( tmp_hypsperc )
    deallocate( tmp_fhru )
    deallocate( tmp_hveg )
    deallocate( tmp_laimax )
end subroutine AWRAL600_main
